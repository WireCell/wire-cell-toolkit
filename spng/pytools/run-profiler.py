import os
import sys
import subprocess
import argparse
import shlex

def main():
    parser = argparse.ArgumentParser(
        description="Run nsys profiler with LD_LIBRARY_PATH wrapper",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with defaults
  python run-profiler.py wire-cell -l stdout config.jsonnet
  
  # With custom nsys arguments
  python run-profiler.py --nsys-args="--trace=cuda,nvtx,osrt --sample=cpu --stats=true --output=my_profile" wire-cell -l stdout config.jsonnet

  # Shorthand
  python run-profiler.py -n "--output=detailed_profile --trace=cuda,nvtx" wire-cell -l stdout config.jsonnet
  
  # Keep the generated script
  python run-profiler.py --keep-script wire-cell -l stdout config.jsonnet
        """
    )
    
    # Single argument for all nsys options
    parser.add_argument(
        "--nsys-args", "-n",
        default="--trace=cuda,nvtx,osrt --sample=cpu --stats=true --output=detailed_profile",
        help='Arguments to pass to nsys profile (default: "--trace=cuda,nvtx,osrt --sample=cpu --stats=true --output=detailed_profile")'
    )
    
    parser.add_argument("--keep-script", action="store_true", help="Keep the generated shell script")
    
    # Capture remaining arguments as the command to run
    parser.add_argument("command", nargs=argparse.REMAINDER, help="Command to profile")
    
    args = parser.parse_args()
    
    # Remove leading '--' if present (from argparse REMAINDER handling)
    if args.command and args.command[0] == '--':
        args.command = args.command[1:]
    
    if not args.command:
        parser.error("No command specified to profile")
    
    # Parse nsys arguments
    try:
        nsys_args = shlex.split(args.nsys_args)
    except ValueError as e:
        parser.error(f"Invalid --nsys-args format: {e}")
    
    # Create a shell script that prepends the LD_LIBRARY_PATH
    script_name = "run_profile.sh"
    
    try:
        with open(script_name, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f'export LD_LIBRARY_PATH="{os.environ.get("LD_LIBRARY_PATH", "")}"\n')
            f.write('"$@"\n')
        
        # Make script executable
        os.chmod(script_name, 0o755)
        
        print(f"Created {script_name}:")
        print("-" * 60)
        with open(script_name, "r") as f:
            print(f.read())
        print("-" * 60)
        
        # Build nsys command
        nsys_cmd = [
            "nsys", "profile",
            *nsys_args,
            f"./{script_name}",
            *args.command
        ]
        
        print(f"\nExecuting: {' '.join(nsys_cmd)}\n")
        print("=" * 60)
        

        # Run the nsys command
        result = subprocess.run(nsys_cmd)
        
        print("=" * 60)
        
        # Extract output filename from nsys_args
        output_name = "detailed_profile"
        for arg in nsys_args:
            if arg.startswith("--output="):
                output_name = arg.split("=", 1)[1]
            elif arg == "--output" and nsys_args.index(arg) + 1 < len(nsys_args):
                output_name = nsys_args[nsys_args.index(arg) + 1]
        
        print(f"\n✓ Profiling complete. Output: {output_name}.nsys-rep")
        
        sys.exit(result.returncode)
        
    except FileNotFoundError:
        print("Error: nsys command not found. Make sure NVIDIA Nsight Systems is installed.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        # Cleanup the script unless --keep-script is specified
        if not args.keep_script and os.path.exists(script_name):
            os.remove(script_name)
            print(f"Removed temporary script: {script_name}")

if __name__ == "__main__":
    main()