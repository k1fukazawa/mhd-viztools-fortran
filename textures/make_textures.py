# Schematic surface images (band structure) of Jupiter and Saturn. Not real images (see README)
import numpy as np
from PIL import Image, ImageFilter
W, H = 1024, 512
rng = np.random.default_rng(7)
lon = (np.arange(W) + 0.5) / W * 2 * np.pi - np.pi          # -pi..pi, longitude 0 at the center
lat = 90 - (np.arange(H) + 0.5) / H * 180                     # top = north pole
LON, LAT = np.meshgrid(lon, lat)

def turb(amp, nmode=24, seed=0):
    r = np.random.default_rng(seed); f = np.zeros_like(LON)
    for k in range(nmode):
        m = r.integers(1, 40); n = r.uniform(0.5, 12)
        f += r.uniform(0.3, 1) / (1 + 0.15 * m) * np.sin(m * LON + n * np.radians(LAT) * 6 + r.uniform(0, 2 * np.pi))
    return amp * f / nmode * 4

def bands(table, latw):
    lat_t = np.array([t[0] for t in table]); col = np.array([t[1] for t in table], float)
    img = np.zeros((H, W, 3))
    for c in range(3):
        img[..., c] = np.interp(-latw, -lat_t, col[:, c])
    return img

# ---- Jupiter (north to south: latitude, RGB)
jt = [(90,(120,112,105)),(60,(140,128,112)),(45,(165,150,128)),(39,(190,170,140)),(35,(160,128,98)),
      (31,(212,196,168)),(27,(150,110,80)),(24,(155,112,82)),(21,(225,212,190)),(18,(222,208,186)),
      (15,(150,100,70)),(10,(140,92,64)),(7,(220,200,170)),(0,(236,222,196)),(-7,(225,205,172)),
      (-9,(160,108,74)),(-17,(148,98,68)),(-20,(228,214,192)),(-25,(230,216,196)),(-28,(165,125,95)),
      (-32,(170,132,100)),(-36,(214,198,172)),(-42,(175,150,125)),(-55,(150,138,120)),(-90,(118,110,104))]
latw = LAT + turb(1.6, seed=1) * (1 + 0.6 * np.cos(np.radians(LAT)))
img = bands(jt, latw)
img *= (1 + 0.06 * turb(1.0, 40, seed=2))[..., None]
# Great Red Spot (22 deg S, longitude 0)
d = ((np.degrees(LON) - 0) / 9.0) ** 2 + ((LAT + 22.5) / 5.5) ** 2
grs = np.exp(-d ** 2)
img = img * (1 - grs[..., None]) + np.array([200, 110, 80]) * grs[..., None]
ring = np.exp(-((np.sqrt(d) - 1.15) / 0.18) ** 2) * 0.35
img = img * (1 - ring[..., None]) + np.array([235, 222, 200]) * ring[..., None]
im = Image.fromarray(np.clip(img, 0, 255).astype(np.uint8)).filter(ImageFilter.GaussianBlur(0.8))
im.save('jupiter.jpg', quality=88)

# ---- Saturn (low-contrast bands)
st = [(90,(150,140,110)),(75,(170,160,125)),(60,(196,180,140)),(45,(214,196,152)),(35,(200,180,138)),
      (25,(222,204,160)),(15,(214,192,146)),(5,(232,214,170)),(0,(236,220,178)),(-5,(230,212,168)),
      (-15,(216,194,150)),(-25,(224,206,162)),(-35,(204,184,142)),(-45,(214,196,152)),(-60,(196,180,140)),
      (-75,(172,160,125)),(-90,(150,140,110))]
latw = LAT + turb(0.8, seed=3)
img = bands(st, latw) * (1 + 0.03 * turb(1.0, 30, seed=4))[..., None]
Image.fromarray(np.clip(img, 0, 255).astype(np.uint8)).filter(ImageFilter.GaussianBlur(1.0)).save('saturn.jpg', quality=88)
