# physics - Double Pendulum

# Project structure

```
📁 physics - Double pendulum
├─ 📄 CPI_ProjetPendule.pdf             # Project informations
├─ 📄 double_pendule_Euler.jl           # Double pendulum simulation with euler
├─ 📄 double_pendule_Euler.mp4
├─ 📄 double_pendule_Runge-Kutta.jl     # Double pendulum simulation with RK4
├─ 📄 double_pendule_Runge-Kutta.mp4
├─ 📄 Manifest.toml
├─ 📄 Project.toml
├─ 📄 README.md
├─ 📁 .vscode
│   └─ 📄 settings.json                 # vscode settings
├─ 📁 analyse
│   ├─ 📄 analyse.trk                   # tracker file
│   ├─ 📄 analyse_video.jl              # Retrieves angles from video in julia
│   ├─ 📄 angles.csv                    # All angles at each frames
│   ├─ 📄 double_pendule.mp4            # Double pendulum video with tracked angles
│   ├─ 📄 find_coordonnee.jl            # Obtain the pivot coordinates in px
│   └─ 📁 video
│       └─ 📄 First_Video_2s.mp4        # Double pendulum real video
└─ 📁 Theorie
    ├─ 📄 graph.jl                      # Schema of the double pendulum
    ├─ 📄 main.tex                      # Theory
    ├─ 📄 Physique_Double_pendule.pdf   # Theory
    └─ 📁 img
        ├─ 📄 Force_m1.png
        ├─ 📄 Force_m2.png
        ├─ 📄 NelderMead.png
        ├─ 📄 RK4.png
        ├─ 📄 shema_double_pendule.jpg
        └─ 📄 shema_double_pendule.svg
```

# Setup

Clone this repository and run
```bash
# Lauch julia
julia
# Open the package manager
]
# Activate the environment
activate .
# Download all dependencies
instantiate
```

# Parameter measurement

- l1 = 91.74 mm
- l2 = 69.33 mm
- Total size of the double pendulum = 20,32 cm
- Size without base = 19 cm

# Final results
## Ecin
https://github.com/ffavre13/Physique_double-pendule/blob/43f762683ff04f6097e784796058eed6520f3c09/results/E_cin.png

## Epot
https://github.com/ffavre13/Physique_double-pendule/blob/43f762683ff04f6097e784796058eed6520f3c09/results/E_pot.png

## Ecin & Epot
https://github.com/ffavre13/Physique_double-pendule/blob/43f762683ff04f6097e784796058eed6520f3c09/results/E_cin-E_pot.png

## Etot
https://github.com/ffavre13/Physique_double-pendule/blob/43f762683ff04f6097e784796058eed6520f3c09/results/E_tot.png

## positions
https://github.com/ffavre13/Physique_double-pendule/blob/43f762683ff04f6097e784796058eed6520f3c09/results/positions.png

## error with NRMSE
- 0.24982410228574903

## Video with the comparison
https://github.com/ffavre13/Physique_double-pendule/blob/43f762683ff04f6097e784796058eed6520f3c09/results/double_pendule_Runge-Kutta_comparison.mp4

## Final video with prediction
https://github.com/ffavre13/Physique_double-pendule/blob/43f762683ff04f6097e784796058eed6520f3c09/results/double_pendule_Runge-Kutta.mp4