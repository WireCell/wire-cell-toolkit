# Mathematical Derivation of Gaussian Weight Calculation

## 1. Basic Setup and Goal

Induced current of drifted charge is computed by convolution of continuous charge distribution and discrete field responses. 
The fundamental problem is to determine how to distribute charge between two adjacent impact wire positions (where the field response is simulated; usually every 1/10th of the wire pitch), in order to effectively apply a linear interpolation of the field responses between the two impact positions, thus maintaining the continuity and accuracy of the simulated waveforms across impact positions along wire pitch orientation for any given charge. 

### Initial Conditions:
- Gaussian distribution centered at μ with standard deviation σ
- Two adjacent wires at positions x₁ and x₂
- Total charge Q in the region between x₁ and x₂

## 2. Mathematical Framework

### The Gaussian Distribution:
```
G(x) = (1/√(2πσ²)) * exp(-(x-μ)²/(2σ²))
```

### Linear Interpolation Position:
For any point x between two wires:
```
w_linear(x) = (x₂ - x)/(x₂ - x₁)  [for wire at x₁]
w_linear(x) = (x - x₁)/(x₂ - x₁)  [for wire at x₂]
```

## 3. Weight Derivation Steps

### Step 1: Consider Charge Conservation
The total charge Q between x₁ and x₂ is:
```
Q = ∫(x₁ to x₂) G(x) dx
```

### Step 2: Center of Charge
The average position of charge (first moment) is:
```
x_avg = (1/Q) * ∫(x₁ to x₂) x*G(x) dx
```

### Step 3: Weight Formula Derivation
The weight formula:
```
w(x) = -σ²/(x₂-x₁) * (G₂-G₁) / Q + (μ-x₂)/(x₁-x₂)
```

This comes from combining:
1. The normalized derivative of the Gaussian (rate of change term)
2. Linear interpolation based on the center position

## 4. Physical Interpretation

The weight calculation balances two physical aspects:
1. The shape of the charge distribution (Gaussian term)
2. The geometric position of the charge relative to the wires (linear term)


```svg
<svg viewBox="0 0 800 900" xmlns="http://www.w3.org/2000/svg">
    <!-- Three panel layout -->
    
    <!-- Panel 1: Basic Gaussian and Wire Setup -->
    <g transform="translate(0,0)">
        <!-- Background -->
        <rect x="50" y="20" width="700" height="250" fill="#f8f8f8" stroke="#ccc"/>
        <text x="400" y="45" text-anchor="middle" font-size="16" font-weight="bold">1. Basic Setup</text>
        
        <!-- Coordinate system -->
        <line x1="100" y1="220" x2="700" y2="220" stroke="black" stroke-width="2"/>
        <line x1="100" y1="220" x2="100" y2="70" stroke="black" stroke-width="2"/>
        
        <!-- Gaussian curve -->
        <path d="M 150 220 
                 C 150 220, 250 210, 300 180 
                 S 350 80, 400 70
                 S 450 80, 500 180
                 S 550 210, 650 220" 
              fill="none" stroke="blue" stroke-width="2"/>
        
        <!-- Wires -->
        <line x1="300" y1="215" x2="300" y2="225" stroke="black" stroke-width="3"/>
        <line x1="500" y1="215" x2="500" y2="225" stroke="black" stroke-width="3"/>
        <text x="300" y="240" text-anchor="middle">x₁</text>
        <text x="500" y="240" text-anchor="middle">x₂</text>
        
        <!-- Center line -->
        <line x1="400" y1="70" x2="400" y2="220" stroke="red" stroke-dasharray="5,5"/>
        <text x="410" y="90">μ</text>
    </g>
    
    <!-- Panel 2: Charge Integration -->
    <g transform="translate(0,300)">
        <!-- Background -->
        <rect x="50" y="20" width="700" height="250" fill="#f8f8f8" stroke="#ccc"/>
        <text x="400" y="45" text-anchor="middle" font-size="16" font-weight="bold">2. Charge Integration</text>
        
        <!-- Coordinate system -->
        <line x1="100" y1="220" x2="700" y2="220" stroke="black" stroke-width="2"/>
        <line x1="100" y1="220" x2="100" y2="70" stroke="black" stroke-width="2"/>
        
        <!-- Gaussian curve -->
        <path d="M 150 220 
                 C 150 220, 250 210, 300 180 
                 S 350 80, 400 70
                 S 450 80, 500 180
                 S 550 210, 650 220" 
              fill="none" stroke="blue" stroke-width="2"/>
        
        <!-- Integration area -->
        <path d="M 300 220 L 300 180 C 350 150, 375 80, 400 70 C 425 80, 450 150, 500 180 L 500 220 Z" 
              fill="rgba(0,255,0,0.2)" stroke="green"/>
        
        <!-- Integration label -->
        <text x="400" y="150" text-anchor="middle">Q = ∫G(x)dx</text>
    </g>
    
    <!-- Panel 3: Weight Calculation -->
    <g transform="translate(0,600)">
        <!-- Background -->
        <rect x="50" y="20" width="700" height="250" fill="#f8f8f8" stroke="#ccc"/>
        <text x="400" y="45" text-anchor="middle" font-size="16" font-weight="bold">3. Weight Components</text>
        
        <!-- Coordinate system -->
        <line x1="100" y1="220" x2="700" y2="220" stroke="black" stroke-width="2"/>
        <line x1="100" y1="220" x2="100" y2="70" stroke="black" stroke-width="2"/>
        
        <!-- Gaussian derivative component -->
        <path d="M 300 180 C 350 150, 375 80, 400 70" fill="none" stroke="purple" stroke-width="2"/>
        <text x="350" y="140" fill="purple">Gaussian derivative term</text>
        
        <!-- Linear interpolation component -->
        <line x1="300" y1="180" x2="500" y2="180" stroke="orange" stroke-width="2"/>
        <text x="400" y="170" fill="orange">Linear position term</text>
        
        <!-- Formula -->
        <text x="400" y="100" text-anchor="middle" font-size="12">
            w(x) = -σ/(x₂-x₁) * (G₂-G₁)/(√2π) / Q + (μ-x₂)/(x₁-x₂)
        </text>
    </g>
</svg>

```

The mathematical interpretation of the weight calculation can be understood in three key parts:

1. **Gaussian Charge Distribution**
   - The fundamental distribution is Gaussian, representing the diffused charge
   - This accounts for the physical nature of charge diffusion in the detector
   - The distribution has a center μ and width σ

2. **Charge Integration**
   - Between any two wires (x₁ and x₂), we have a total charge Q
   - This charge must be conserved in our weighting scheme
   - The integral of the Gaussian between the wires gives us Q:
     ```
     Q = ∫(x₁ to x₂) G(x) dx
     ```

3. **Weight Calculation Components**
   The weight formula combines two essential physical aspects:

   a) **Gaussian Derivative Term**: `-σ/(x₂-x₁) * (G₂-G₁)/(√2π) / Q`
      - Represents how quickly the charge distribution changes
      - Accounts for the shape of the Gaussian
      - Normalized by total charge Q to maintain conservation

   b) **Linear Position Term**: `(μ-x₂)/(x₁-x₂)`
      - Provides basic geometric interpolation
      - Ensures smooth transition between wires
      - Based on the center position relative to wire positions

The key insights behind this formulation are:

1. **Physical Motivation**:
   - The charge distribution is continuous but must be measured at discrete points
   - The weighting should reflect both the distribution shape and position
   - Total charge must be conserved

2. **Mathematical Properties**:
   - Smooth transition between wires
   - Proper handling of both narrow and wide distributions
   - Conservation of total charge
   - Correct limiting behavior for point-like deposits

3. **Practical Implementation**:
   - Computationally efficient
   - Numerically stable
   - Handles edge cases appropriately

The diagrams show:
1. The basic setup with the Gaussian distribution and wire positions
2. The charge integration between wires
3. The components of the weight calculation

