## This describes how to set up a new ipython kernel on jupyterlab

1. Open a terminal
2. Using your favorite virtual environment manager, make a new venv, then activate it
 
  a. For this demo, I assume we have access to uv
<pre>
  uv venv my_new_venv
  source my_new_venv/bin/activate
</pre>

3. Now, we have to install ipykernel along with any other packages we want

<pre>
  uv pip install ipykernel {other packages}
</pre>

4. We can then install the kernel

<pre>
  uv run ipython kernel install --user --env VENV ${PWD}/my_new_venv/ --name={kernel_name}
</pre>

5. Sometimes your kernel will show up automatically in the 'Launcher' page. If not, log out and back in and it should be there.

At this point you can start using packages installed in the virtual environment we made. You can also continue to add more packages.
If you are in an environment and install another package, you will at least need to _restart_ your kernel. If you can't import the 
package even after restarting the kernel, try logging out and back in and it should work.
