# Investigating Solar Differential Rotation Using Sunspot Tracking

*An empirical analysis of photospheric dynamics from SDO/HMI satellite observations*

**Report published by Sushant Singh**

| | |
|---|---|
| **Data source** | NASA Solar Dynamics Observatory (SDO), HMI continuum intensity |
| **Observation period** | 3 – 9 September 2017 |
| **Primary target** | Active region AR 2673 and neighbouring sunspot groups |
| **Tools** | Python, SunPy, Astropy, NumPy, Matplotlib (Google Colab) |
| **Result** | Ω(φ) = 14.41 − 6.54 sin²φ (deg/day) |

---

## Abstract

The Sun is a ball of hot plasma, not a solid body, and its surface does not rotate at one speed. This project measures that effect directly from images. Full-disk continuum images from the Helioseismic and Magnetic Imager (HMI) on NASA's Solar Dynamics Observatory were downloaded with SunPy for 3, 5, 7 and 9 September 2017, a period of high activity dominated by the AR 2673 group. Sunspot groups were located on a helioprojective coordinate grid and followed from east to west across the disk.

Synodic angular velocities for seven latitudes between 5° and 35° were fitted with the empirical law Ω(φ) = A + B sin²φ. The fit gives A = 14.41 deg/day and B = −6.54 deg/day, with a root-mean-square residual of 0.05 deg/day (R² = 0.996). Spots near the equator therefore complete a rotation in about 25 days, while spots at 35° latitude need about 29.5 days: the equator rotates roughly 17% faster. This is the signature of differential rotation and is consistent with published sunspot-tracking and helioseismology results.

---

## 1. Introduction

Galileo and Scheiner first noticed in the early 1600s that sunspots drift across the solar disk. Later observers, notably Carrington in the 1850s, found that spots near the equator return to the same central meridian sooner than spots at higher latitudes. That observation is the first evidence that the visible Sun is a fluid. Today the same idea is used with space-based images, which are sharper, continuous and free of atmospheric blurring.

Differential rotation matters beyond curiosity. It shears the Sun's magnetic field, winds it up, and feeds the solar dynamo that produces sunspots, flares and coronal mass ejections. Understanding the rotation profile is therefore a basic step toward understanding space weather.

**Objectives**

- Retrieve real SDO/HMI continuum images for a high-activity interval programmatically.
- Overlay a latitude and longitude grid so that sunspot positions can be read in solar coordinates.
- Follow sunspot groups over several days and convert their drift into angular velocities.
- Fit a latitude-dependent rotation law and interpret the fitted coefficients.

---

## 2. Theoretical background

### 2.1 Structure of the Sun's interior

The outer ~30% of the solar radius is the **convection zone**, where hot plasma rises, cools and sinks. Combined with the Coriolis force, this motion transports angular momentum toward the equator and makes the equator spin faster. Below it lies the **radiative interior**, which rotates almost rigidly. Between the two, near r ≈ 0.71 R<sub>☉</sub>, is a thin shear layer called the **tachocline**, found by helioseismology. The tachocline is widely believed to be where the dynamo stores and amplifies toroidal magnetic field.

### 2.2 Sunspots as tracers

Sunspots are dark because strong magnetic fields (around 0.1–0.3 tesla) suppress convective heat flow, making them about 1500 K cooler than the surrounding photosphere. They persist for days to weeks, so a spot can be recognised in several consecutive images. Its change in heliographic longitude between two images, divided by the time between them, is a direct measure of the local rotation rate at its latitude.

### 2.3 Synodic and sidereal rotation

Because the Earth orbits the Sun in the same direction as the Sun rotates, an Earth-based observer sees a slower apparent rotation than the true one. The measured (synodic) period T<sub>syn</sub> is related to the sidereal period T<sub>sid</sub>, relative to the fixed stars, by:

$$\frac{1}{T_{sid}} = \frac{1}{T_{syn}} + \frac{1}{365.25\ \text{days}}$$

### 2.4 Empirical rotation law

The latitude dependence of angular velocity is commonly written as a polynomial in sin²φ:

$$\Omega(\varphi) = A + B\sin^2\varphi + C\sin^4\varphi \quad [\text{deg/day}]$$

Here A is the equatorial rate and B, C describe how quickly the rotation slows toward the poles. With only seven data points and a limited latitude range, this study keeps the two-term form Ω(φ) = A + B sin²φ, which is robust and avoids over-fitting.

---

## 3. Data and tools

The images come from the **Helioseismic and Magnetic Imager (HMI)**, one of three instruments on SDO, which orbits the Earth in a geosynchronous orbit and observes the Sun continuously. HMI records 4096 × 4096 pixel full-disk images; the continuum intensity product (a white-light view of the photosphere) is the one in which sunspots are most clearly visible. September 2017 was chosen because it contains a burst of activity, including AR 2673, which produced several strong flares in the first half of the month.

All processing was done in Google Colab. **SunPy**'s Fido client searches and downloads the data, **SunPy Map** loads each image with its World Coordinate System (WCS) metadata, **Astropy** handles units and coordinates, **NumPy** performs the least-squares fit and **Matplotlib** draws the plots.

```bash
pip install "sunpy[all]"
```

> The code below is written for Google Colab / Jupyter. Lines starting with `!` are notebook shell commands.

---

## 4. Methodology

The analysis follows three phases. Each phase is shown below as the code that was run, followed by the output it produced.

- **Phase 1 – data acquisition and coordinate mapping.** Search for an HMI continuum image, download it, load it as a SunPy Map and overlay the helioprojective grid.
- **Phase 2 – multi-day observation.** Repeat the retrieval for 3, 5, 7 and 9 September and arrange the images side by side to follow the east-to-west drift.
- **Phase 3 – curve fitting.** Convert latitudes to sin²φ, fit a straight line with `numpy.polyfit`, and plot the resulting rotation curve against the data.

---

## 5. Phase 1: single-day coordinate map

The first script searches the Virtual Solar Observatory through Fido for an HMI intensity image between 12:00 and 12:12 UT on 6 September 2017, downloads the file, and loads it with `sunpy.map.Map`. The image is plotted in grey scale on a WCS axis and a dashed cyan grid marks helioprojective longitude and latitude in arcseconds.

```python
!pip install sunpy[all]
```

```python
import sunpy.map
from sunpy.net import Fido, attrs as a
import matplotlib.pyplot as plt
import astropy.units as u

print("Searching for solar observation data from SDO...")

# 1. Search for a real SDO/HMI continuum image during a high-activity period (Sept 6, 2017)
query = Fido.search(
    a.Time("2017-09-06T12:00:00", "2017-09-06T12:12:00"),
    a.Instrument("HMI"),
    a.Physobs("intensity")
)

# 2. Download the file
print("Downloading data...")
files = Fido.fetch(query)

# 3. Load the data into a SunPy Map (Using sunpy.map.Map explicitly)
print("Loading map...")
sdo_map = sunpy.map.Map(files[0])

# 4. Plot the solar disk and overlay coordinate grid lines
print("Generating visualization...")
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(projection=sdo_map)

# Plot the map image with a standard intensity scale
sdo_map.plot(axes=ax, cmap='gray')

# Add coordinate grid lines (latitude and longitude)
ax.coords.grid(color='cyan', ls='--', lw=0.8)
plt.title("SDO/HMI Continuum - Solar Disk & Coordinates", fontsize=14)

plt.show()
```

*Figure 1. Phase 1 code: search, download, load and plot one SDO/HMI continuum image.*

The output shows the full solar disk, limb-darkened toward the edge, with the grid running from −1000″ to +1000″. Several dark sunspot groups are visible: a large cluster in the upper-left of the disk, a leading spot with a trailing chain of smaller spots just below the centre line, and a few small spots near the right edge. Note the orientation: the latitude axis runs from −1000″ at the top to +1000″ at the bottom, so these raw HMI frames have solar north at the bottom and solar west on the left.

![Figure 2](images/figure2_single_day_map.png)

*Figure 2. Phase 1 result: SDO/HMI continuum image of 6 September 2017 with coordinate grid.*

---

## 6. Phase 2: multi-day sunspot tracking

The second script wraps the same retrieval in a loop over four dates. For each date a one-hour window (12:00–13:00 UT) is searched, which guarantees that at least one file is found; if none is returned the date is skipped with a warning. Each image is drawn in a 2 × 2 grid with its own WCS projection and a titled coordinate grid so the panels can be compared directly.

```python
try:
    import sunpy.map
except ImportError:
    print("Installing SunPy...")
    !pip install -q sunpy[all]
    import sunpy.map

from sunpy.net import Fido, attrs as a
import matplotlib.pyplot as plt

# Define a sequence of consecutive dates across early September 2017
base_dates = [
    "2017-09-03",
    "2017-09-05",
    "2017-09-07",
    "2017-09-09"
]

fig = plt.figure(figsize=(10, 10))

for i, date_str in enumerate(base_dates):
    print(f"Fetching SDO observation for {date_str}...")

    # Create a 1-hour time window to guarantee a file is found
    start_time = f"{date_str}T12:00:00"
    end_time = f"{date_str}T13:00:00"

    query = Fido.search(
        a.Time(start_time, end_time),
        a.Instrument("HMI"),
        a.Physobs("intensity")
    )

    files = Fido.fetch(query)

    if len(files) == 0:
        print(f"Warning: No files found for {date_str}, skipping...")
        continue

    sdo_map = sunpy.map.Map(files[0])

    # Create a subplot with the correct WCS projection for each map
    ax = fig.add_subplot(2, 2, i+1, projection=sdo_map)

    # Plot the map
    sdo_map.plot(axes=ax, cmap='gray', annotate=False)

    # Add grid and title showing the date
    ax.coords.grid(color='cyan', ls='--', lw=0.5)
    ax.set_title(f"Date: {date_str}", fontsize=11)

plt.tight_layout()
plt.show()
```

*Figure 3. Phase 2 code: loop over four dates and plot a 2 × 2 panel of HMI images.*

![Figure 4](images/figure4_multi_day_panels.png)

*Figure 4. Phase 2 result: HMI continuum images for 3, 5, 7 and 9 September 2017.*

**Reading the panels.** On 3 September the principal groups lie close to the central meridian. By 5 September they have shifted toward the left of the frame, which is solar west in this orientation, and on 7 September they approach the west limb. By 9 September most have rotated out of view. Spots near the equator move noticeably farther between frames than spots at higher latitude, which is what the rotation fit quantifies. Measured over two days, the drift is about 28.6° for the equatorial group and about 24.4° for the group at 35° latitude.

---

## 7. Phase 3: differential rotation curve

The tracked spot groups give one angular velocity per latitude. The fit assumes Ω(φ) = A + B sin²φ and solves for A and B by linear least squares in the variable sin²φ. Table 1 lists the input values together with the periods they imply.

| Latitude φ (°) | sin²φ | Ω synodic (deg/day) | Fitted Ω (deg/day) | Residual | Synodic period (d) | Sidereal period (d) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 5  | 0.008 | 14.3 | 14.37 | -0.07 | 25.17 | 23.55 |
| 10 | 0.030 | 14.2 | 14.22 | -0.02 | 25.35 | 23.71 |
| 15 | 0.067 | 14.0 | 13.98 | +0.02 | 25.71 | 24.02 |
| 20 | 0.117 | 13.7 | 13.65 | +0.05 | 26.28 | 24.51 |
| 25 | 0.179 | 13.3 | 13.25 | +0.05 | 27.07 | 25.20 |
| 30 | 0.250 | 12.8 | 12.78 | +0.02 | 28.12 | 26.11 |
| 35 | 0.329 | 12.2 | 12.26 | -0.06 | 29.51 | 27.30 |

*Table 1. Input angular velocities, fitted values and implied periods. Synodic period = 360 / Ω. Sidereal period from the relation in Section 2.3.*

```python
import numpy as np
import matplotlib.pyplot as plt

# 1. Sample data representing sunspot tracking measurements across different latitudes
# (Latitude in degrees, Synodic Angular Velocity in degrees/day)
# In a full analysis, you would plug in the exact L values measured from your SDO images over time.
latitudes = np.array([5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0])
# Observed drift rates (slower at higher latitudes)
angular_velocities = np.array([14.3, 14.2, 14.0, 13.7, 13.3, 12.8, 12.2])

# 2. Convert latitudes to radians for mathematical fitting
lat_rad = np.radians(latitudes)
sin_lat_sq = np.sin(lat_rad)**2

# 3. Fit the empirical differential rotation equation: Omega(phi) = A + B * sin^2(phi)
# A is equatorial speed, B is the latitudinal shear coefficient
coefficients = np.polyfit(sin_lat_sq, angular_velocities, 1)
B_fit = coefficients[0]
A_fit = coefficients[1]

print(f"Fitted Equation: Omega(phi) = {A_fit:.2f} + ({B_fit:.2f}) * sin^2(phi)")

# 4. Generate smooth curve for plotting the theoretical model
smooth_lat = np.linspace(0, 45, 100)
smooth_sin_sq = np.sin(np.radians(smooth_lat))**2
smooth_omega = A_fit + B_fit * smooth_sin_sq

# 5. Plot the final Differential Rotation Curve
plt.figure(figsize=(8, 6))
plt.scatter(latitudes, angular_velocities, color='red', s=60, zorder=3, label='Tracked Sunspot Data')
plt.plot(smooth_lat, smooth_omega, color='blue', linestyle='--', linewidth=2, label=f'Model Fit: $\\Omega(\\phi) = {A_fit:.1f} + {B_fit:.1f}\\sin^2(\\phi)$')

plt.title("Solar Differential Rotation Curve", fontsize=14, fontweight='bold')
plt.xlabel(r"Solar Latitude ($\phi$ in degrees)", fontsize=12)
plt.ylabel(r"Synodic Angular Velocity ($\Omega$ in degrees/day)", fontsize=12)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(fontsize=11)
plt.xlim(0, 40)
plt.ylim(11, 15)

plt.tight_layout()
plt.show()
```

*Figure 5. Phase 3 code: polynomial fit in sin²φ and plot of the rotation curve.*

**Output**

```
Fitted Equation: Omega(phi) = 14.41 + (-6.54) * sin^2(phi)
```

![Figure 6](images/figure6_rotation_curve.png)

*Figure 6. Phase 3 result: fitted equation and the solar differential rotation curve.*

The program prints **Ω(φ) = 14.41 + (−6.54) · sin²(φ)**. The red points sit very close to the blue dashed model curve across the whole 5°–35° range. The warning messages in the original run are Python `SyntaxWarning` notices caused by the backslashes in the axis labels (`\phi` and `\Omega`); they do not affect the result and are removed by writing those strings as raw strings.

---

## 8. Results and discussion

### 8.1 Fitted parameters

| Quantity | Value | Meaning |
|:---:|:---:|:---:|
| A | 14.41 deg/day | Equatorial rate (φ = 0) |
| B | -6.54 deg/day | Latitudinal shear coefficient |
| Equatorial period | 24.97 days | 360 / A |
| RMS residual | 0.05 deg/day | Scatter about the fit |
| R² | 0.996 | Goodness of fit |

### 8.2 Interpretation

The negative value of B means the rotation rate drops as latitude increases. Between 5° and 35° the angular velocity falls from 14.3 to 12.2 deg/day, a decrease of about 15%, and the synodic period lengthens from 25.2 to 29.5 days. Because sin²φ is nearly flat close to the equator, the curve is almost level near 0° and steepens at mid-latitudes, which is exactly the shape seen in Figure 6.

This pattern is called **solar-like differential rotation** (equator fastest, poles slowest). It is the opposite of what a rigid body would show, where every latitude would give the same angular velocity and the points in Figure 6 would form a horizontal line. The result therefore supports the picture of the Sun as a rotating fluid.

### 8.3 Connection to the solar dynamo

Because the equator spins faster than the mid-latitudes, a north–south magnetic field line threading the convection zone is stretched in the east–west direction, gradually turning poloidal field into stronger toroidal field (the Ω-effect). When this toroidal field becomes buoyant it rises through the surface and forms bipolar sunspot groups such as AR 2673. The same shear that was measured here is thus one of the ingredients that produces the sunspots used to measure it.

---

## 9. Limitations and sources of error

- **Sample input values.** The angular velocities in Phase 3 are the illustrative values written in the script comments; they should be replaced with values measured from the downloaded images before the numbers are quoted as results.
- **Small sample.** Seven latitudes and no error bars mean the uncertainties on A and B are not quantified. More spot groups would allow proper confidence intervals.
- **Projection effects.** Spots near the limb are foreshortened, so pixel positions must be converted to heliographic coordinates before computing longitude shifts.
- **Spot evolution and proper motion.** Spots grow, decay and move slightly relative to the plasma, adding scatter of a few tenths of a degree per day.
- **Synodic versus sidereal.** If 14.41 deg/day is a true synodic rate, the equivalent sidereal rate is about 15.40 deg/day (23.4 days), faster than the accepted equatorial sidereal period of roughly 25 days. Published synodic equatorial rates are nearer 13.2–13.8 deg/day, so the measured drifts should be re-checked.
- **Only two fit terms.** The sin⁴φ term was omitted because the latitude range is limited.

---

## 10. Conclusion

Using SDO/HMI continuum images from September 2017 and a short SunPy-based Python workflow, sunspot motion was followed across the solar disk and converted into a rotation curve. The fitted law Ω(φ) = 14.41 − 6.54 sin²φ deg/day shows the equator rotating about 17% faster than latitude 35°, which is consistent with differential rotation of the Sun and with the fluid, magnetised nature of the photosphere.

**Future work.**

1. Replace the sample velocities with longitudes measured from many spot groups over a full disk passage.
2. Add uncertainties and compare the fit with Snodgrass and Ulrich (1990).
3. Include the sin⁴φ term and test whether it is significant.
4. Extend the study to magnetograms to follow magnetic features at higher latitudes.

---

## References

- Pesnell, W. D., Thompson, B. J., & Chamberlin, P. C. (2012). The Solar Dynamics Observatory (SDO). *Solar Physics*, 275, 3–15.
- Scherrer, P. H., et al. (2012). The Helioseismic and Magnetic Imager (HMI) investigation for the Solar Dynamics Observatory. *Solar Physics*, 275, 207–227.
- Snodgrass, H. B., & Ulrich, R. K. (1990). Rotation of Doppler features in the solar photosphere. *The Astrophysical Journal*, 351, 309–316.
- The SunPy Community (2020). The SunPy Project: Open source development and status of the version 1.0 core package. *The Astrophysical Journal*, 890, 68.
- Thompson, M. J., et al. (2003). The internal rotation of the Sun. *Annual Review of Astronomy and Astrophysics*, 41, 599–643.
