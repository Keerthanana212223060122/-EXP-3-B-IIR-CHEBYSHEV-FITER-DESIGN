# EXP 3 : IIR-CHEBYSHEV-FITER-DESIGN

# REG NO : 212223060122

## AIM: 
To design an IIR Chebyshev filter  using SCILAB. 

## APPARATUS REQUIRED: 
PC installed with SCILAB. 

## PROGRAM (LPF): 
```python
clc;
clear;
close();

// --- User Inputs ---
wp = input("Enter the pass band frequency (Radians )= ");
ws = input("Enter the stop band frequency (Radians )= ");
alphap = input("Enter the pass band attenuation (dB)= ");
alphas = input("Enter the stop band attenuation(dB)= ");
T = input("Enter the Value of sampling Time= ");

// -------------------------
// Pre-warp (bilinear)
// -------------------------
omegap = (2/T) * tan(wp/2);
mprintf("\nomegap = %f\n", omegap);

omegas = (2/T) * tan(ws/2);
mprintf("omegas = %f\n", omegas);

// -------------------------
// Chebyshev order (analytical)
// -------------------------
N = acosh(sqrt(((10^(0.1*alphas))-1)/((10^(0.1*alphap))-1))) / acosh(omegas/omegap);
mprintf("N = %f\n", N);

N = ceil(N);
mprintf("Round off value of N = %d\n", N);

// -------------------------
// Cutoff & epsilon
// -------------------------
omegac = omegap / (((10^(0.1*alphap)) - 1)^(1/(2*N)));
mprintf("omegac = %f\n", omegac);

Epsilon = sqrt((10^(0.1*alphap)) - 1);
mprintf("Epsilon = %f\n\n", Epsilon);

// -------------------------
// Poles & gain from zpch1
// -------------------------
[pols, gn] = zpch1(N, Epsilon, omegap);
mprintf("Gain = %f\n\n", gn);

// Print poles
mprintf("Poles\n");
for i = 1:length(pols)
    p = pols(i);
    if imag(p) >= 0 then
        mprintf(" - %g + %gi\n", real(p), imag(p));
    else
        mprintf(" - %g - %gi\n", real(p), abs(imag(p)));
    end
end
mprintf("\n");

// -------------------------
// Build analog transfer function H(s)
// -------------------------
s = poly(0, 's'); // define s
hs = poly(gn, 's', 'coeff') / real(poly(pols, 's')); // H(s) = gn / product(s - p_i)
mprintf("Analog Low-Pass Chebyshev Transfer function H(s) =\n");
disp(hs);

// -------------------------
// Bilinear transform -> Digital H(z)
// -------------------------
z = poly(0, 'z'); // define z
Hz = horner(hs, (2/T) * ((z - 1) / (z + 1))); // s -> z
mprintf("\nDigital LPF Transfer function H(z) =\n");
disp(Hz);

// -------------------------
// Frequency response
// -------------------------
HW = frmag(Hz, 512); // frequency samples
w = 0:%pi/511:%pi;

// -------------------------
// Proper colored plot using plot2d
// -------------------------
figure();
plot2d(w/%pi, abs(HW), style=2); // style=2 -> blue solid line
xtitle("Frequency Response of Chebyshev IIR Low-Pass Filter", "Normalized Digital Frequency w", "Magnitude");
gcf().background = color("white"); // set figure background to white
xgrid();
```

## PROGRAM (HPF): 
```python
clc;
close;

// ============================
// User Inputs
// ============================
wp = input('Enter the pass band frequency (Radians )= ');  // Passband edge > Stopband edge
ws = input('Enter the stop band frequency (Radians )= ');
alphap = input('Enter the pass band attenuation (dB)= ');
alphas = input('Enter the stop band attenuation (dB)= ');
T = input('Enter the Value of sampling Time= ');

// ============================
// Pre-warping (Bilinear Transformation)
// ============================
omegap = (2/T)*tan(wp/2);   // Passband edge (analog)
disp(omegap,'omegap=');
omegas = (2/T)*tan(ws/2);   // Stopband edge (analog)
disp(omegas,'omegas=');

// ============================
// Order of the HPF
// ============================
// For HPF: use acosh(omegap/omegas)
N = acosh(sqrt(((10^(0.1*alphas))-1)/((10^(0.1*alphap))-1))) / acosh(omegap/omegas);
disp(N,'N=');
N = ceil(N);
disp(N,'Round off value of N=');

// ============================
// Cutoff frequency
// ============================
omegac = omegap / (((10^(0.1*alphap))-1)^(1/(2*N)));
disp(omegac,'omegac=');

// ============================
// Chebyshev Prototype (Normalized LPF)
// ============================
Epsilon = sqrt((10^(0.1*alphap))-1);
disp(Epsilon,'Epsilon=');

[pols, gn] = zpch1(N, Epsilon, 1);   // normalized LPF prototype at 1 rad/s
disp(gn,'Gain');
disp(pols,'Poles');

s = poly(0,'s');   // Laplace variable
hs = poly(gn,'s','coeff') / real(poly(pols,'s'));
disp(hs,'Analog Normalized Chebyshev LPF');

// ============================
// LPF → HPF Transformation (s -> omegac/s)
// ============================
sh = horner(hs, omegac/s);
disp(sh,'Analog Chebyshev High-Pass Filter');

// ============================
// Bilinear Transformation: s -> (2/T)*((z-1)/(z+1))
// ============================
z = poly(0,'z');   // z-domain variable
Hz = horner(sh, (2/T) * ((z-1)/(z+1)));
disp(Hz,'Digital HPF Transfer function H(Z)=');

// ============================
// Frequency Response
// ============================
HW = frmag(Hz, 512);
w = 0:%pi/511:%pi;

plot(w/%pi, abs(HW));
xlabel('Normalized Digital Frequency ×π rad/sample');
ylabel('Magnitude');
title('Frequency Response of Chebyshev IIR High-Pass Filter');
xgrid();
```

## OUTPUT (LPF) : 

<img width="1919" height="988" alt="Screenshot 2025-10-22 121322" src="https://github.com/user-attachments/assets/8f1a8768-9445-4d7f-9e28-2e4b398f76a1" />

<img width="1679" height="844" alt="Screenshot 2025-10-22 121347" src="https://github.com/user-attachments/assets/f61d8dd4-b746-41b2-b2d0-32d0eb6a28ad" />

## OUTPUT (HPF) : 

<img width="1919" height="992" alt="Screenshot 2025-10-22 122207" src="https://github.com/user-attachments/assets/66639bdb-b71e-4524-99e3-87f97d687412" />

<img width="1664" height="828" alt="Screenshot 2025-10-22 121640" src="https://github.com/user-attachments/assets/1b9fb2b6-13cc-4729-be13-8debb03f80b9" />

## RESULT: 

IIR-CHEBYSHEV-FITER-DESIGN OF LPF AND HPF IS SUCCESSFULLY EXECUTED.
