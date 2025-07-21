# 🌈 100 Gbps WDM Optical Link using Lumerical

This project simulates a **10-channel, 100 Gbps Wavelength Division Multiplexing (WDM) optical communication system** using **Lumerical INTERCONNECT**, and analyzes system performance using **MATLAB eye diagram analysis**.

---

## ⚙️ Project Overview

- **Modulation**: 10 × 10 Gbps NRZ using Mach-Zehnder Modulators (MZMs)
- **Laser Sources**: CW lasers starting from **193.6 THz**
- **Multiplexing**: 10×1 Optical Multiplexer (WDM)
- **Medium**: Optical Waveguide
- **Receiver**: Photodiode
- **Analysis**: Eye Diagram analysis (Q-factor, BER) using MATLAB

---

## 🎯 Objectives

- Demonstrate high-bandwidth transmission using WDM
- Improve spectral efficiency via channel multiplexing
- Evaluate quality of signal using Q-factor and BER from eye diagrams

---

## 📁 Folder Structure

| Folder       | Description                           |
|--------------|----------------------------------------|
| `circuit/`   | Lumerical design files (`.lsf`, `.icp`) |
| `analysis/`  | MATLAB script for eye diagram analysis |
| `results/`   | Eye diagram screenshots, WDM plots     |
| `docs/`      | (Optional) Project report or summary   |

---




## 📈 Sample Output

Here are the performance metrics comparing the **Photonic Circuit** (simulated in Lumerical) with a **Traditional Electrical Circuit**, both using 10 Gbps NRZ signaling:

| **Metric**                              | **Photonic Circuit** | **Traditional Circuit** |
|----------------------------------------|-----------------------|--------------------------|
| Signal Type                             | 10 Gbps NRZ           | 10 Gbps NRZ              |
| Bit Period (ps)                         | 100                   | 100                      |
| Eye Height (mA)                         | 0.22                  | *0.15                    |
| Eye Amplitude (mA)                      | 2.44                  | *2.0                     |
| Extinction Ratio (dB)                  | 20.23                 | *10                      |
| Approximated Q-Factor                   | 26.81                 | *6                       |
| Signal-to-Noise Ratio (dB)              | 31.55                 | *20                      |
| Eye Threshold for Crossings (mA)        | 1.24                  | *1.0                     |
| Total Number of Crossings Detected      | 3709                  | *3500                    |
| Min Crossing Time (ps)                  | -24.38                | *-25                     |
| Max Crossing Time (ps)                  | 50                    | *50                      |
| Mean Crossing Time (ps)                 | 14.1                  | *15                      |
| Peak-to-Peak Jitter (ps)                | 54.69                 | *60                      |
| RMS Jitter (ps)                         | 17.4                  | *20                      |
| Latency (ps)                            | 215                   | 550                      |

> * indicates estimated or typical values for conventional electrical circuits from literature or prior simulation studies.

---

## 🚀 How to Run

1. Open the `.lsf` or `.icp` file in **Lumerical INTERCONNECT**.
2. Run the simulation for all 10 channels.
3. Export eye diagram data for each channel.
4. Use `eye_diagram_analysis.m` in MATLAB to compute Q-factor and BER.

---

## 🧠 Tech Stack

- [ Ansys Lumerical INTERCONNECT](https://www.ansys.com/en-in/products/optics/interconnect)
- [MATLAB](https://www.mathworks.com/products/matlab.html)

---

