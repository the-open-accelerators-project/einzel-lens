r"""
================================================================================
 EINZELLINSE
==============================================================================

requirements: pip install numpy matplotlib

--------------------------------------------------------------------------------
WAS DIESE RECHNUNG NICHT KANN
--------------------------------------------------------------------------------
* Elektronen weit außen (mehr als ca. 1/3 des Radius): Das kann nur die große (Laplace) Simulation.
* Energieunschärfe: 24 +- 0,5 eV verschiebt den Fokus bei -16 V um ±2 mm.
* !Erdmagnetfeld!: lenkt 24-eV-Elektronen auf 33 mm schon ~1,5 mm ab.
==================
"""
import numpy as np
import matplotlib.pyplot as plt

# =========================== EINGABEN ============================
D = 10.0        # Bohrung = Innendurchmesser der Ringe
L = 10.0        # Länge des MITTLEREN Rings (die äußeren sind hier egal,
                # solange sie mindestens ~D lang sind)
g = 2.0         # Spaltbreite zwischen den Ringen
E = 24.0        # Strahlenergie in eV (= Kathodenspannung, nur ohne Minus)
U_BILD = -16.0  # für diese Spannung wird ein Bild gezeichnet
SPANNUNGEN = [
  -50, -49, -48, -47, -46, -45, -44, -43, -42, -41,
  -40, -39, -38, -37, -36, -35, -34, -33, -32, -31,
  -30, -29, -28, -27, -26, -25, -24, -23, -22, -21,
  -20, -19, -18, -17, -16, -15, -14, -13, -12, -11,
  -10, -9, -8, -7, -6 ,-5, -4, -3, -2, -1, -0.5,  0.00001,  5, 10, 15, 20, 30, 45, 60
]
# =============================================================================

a = D / 2                       # Rohrradius
w = 1.318                       # Formfaktor der tanh-Näherung (Abschnitt 5)
c = (L + g) / 2                 # Abstand Linsenmitte -> Spaltmitte
dz = 0.01                       # Schrittweite in mm (Abschnitt 7)
z = np.arange(-60, 250, dz)     # alle Stellen z, an denen wir rechnen.
                                # Linsenmitte bei z = 0, Start 60 mm davor.


def S(x):
    """Ein S-förmiger Sprung von -1 auf +1 über einen Spalt (Abschnitt 5)."""
    return a / (w * g) * np.log(np.cosh(w * (x + g/2) / a) / np.cosh(w * (x - g/2) / a))


def bahn(U, r0=1.0):
    """Verfolgt ein Elektron, das im Abstand r0 parallel zur Achse startet.
    Gibt das Potential V(z), die Bahn r(z), den Brennpunkt F, die
    Brennweite f und die Hauptebene H zurück."""
    V = U / 2 * (S(z + c) - S(z - c))      # Potential auf der Achse
    K = V + E                              # Bewegungsenergie in Volt
    if K.min() <= 0:                       # Elektron bleibt stehen -> Spiegel
        return V, None, None, None, None

    # Steigung und Krümmung von K aus den Nachbarwerten
    # Ganz am Rand bleibt 0 stehen, dort ist eh kein Feld mehr.
    K1 = np.zeros_like(K)
    K2 = np.zeros_like(K)
    K1[1:-1] = (K[2:] - K[:-2]) / (2 * dz)
    K2[1:-1] = (K[2:] - 2 * K[1:-1] + K[:-2]) / dz**2

    r, r1 = r0, 0.0                        # Startabstand, Startschräge
    rr = np.empty_like(z)                  # hier merken wir uns die Bahn
    for i in range(len(z)):                # das Daumenkino (Abschnitt 7)
        rr[i] = r
        r2 = -K1[i] / (2 * K[i]) * r1 - K2[i] / (4 * K[i]) * r
        r1 += r2 * dz                      # Schräge anpassen
        r += r1 * dz                       # Position anpassen
        if z[i] >= 60 and i > 0 and z[i - 1] < 60:
            r60, r1_60 = r, r1             # Werte weit hinter der Linse merken

    # Brennpunkt F: erste Stelle hinter der Mitte, an der r von + nach - geht
    kreuz = np.where((z[1:] > 0) & (rr[:-1] > 0) & (rr[1:] <= 0))[0]
    F = None
    if len(kreuz):
        k = kreuz[0]                       # zwischen zwei Punkten linear
        F = z[k] + rr[k] / (rr[k] - rr[k + 1]) * dz   # interpolieren
    f = -r0 / r1_60                        # Brennweite aus der Endschräge
    H = (60 - r60 / r1_60) - f             # Hauptebene = F(Gerade) - f
    return V, rr, F, f, H


# ------------------------- OUTPUT Tabelle für alle Spannungen -----------------------
print(f"Einzellinse: Bohrung {D} mm, Mittelring {L} mm, Spalt {g} mm, Strahl {E} eV")
print("  U (V) |  U/E  | Brennpunkt F | Brennweite f | Hauptebene H   (mm, ab Linsenmitte)")
for U in SPANNUNGEN:
    V, rr, F, f, H = bahn(U)
    if rr is None:
        print(f"{U:7.0f} | {U/E:+.2f} | Elektron wird reflektiert (|U| >= E)")
    else:
        Ftxt = f"{F:8.1f}" if F else "  > 250 "
        print(f"{U:7.0f} | {U/E:+.2f} |   {Ftxt}   |   {f:8.1f}   |   {H:+6.1f}")

# ------------------------------- Bild zeichnen; Dieser Abschnitt ist mit Hilfe von Claude geschrieben -------------------------------
V, rr, F, f, H = bahn(U_BILD)
fig, (o, u) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
o.plot(z, V + E, color="#4a3aa7")
o.axhline(E, color="grey", lw=0.5, ls="--")
o.set_ylabel("K = Bewegungsenergie (eV)")
o.set_title(f"U = {U_BILD:.0f} V am Mittelring, Strahl {E:.0f} eV")
for s in (-1, 1):                                   # Spalte markieren
    o.axvspan(s * c - g/2, s * c + g/2, color="#e34948", alpha=0.15)
    u.axvspan(s * c - g/2, s * c + g/2, color="#e34948", alpha=0.15)
u.plot(z, rr, color="#e34948", label="Bahn (Start r = 1 mm)")
u.axhline(0, color="black", lw=0.6)
if F:
    u.axvline(F, color="grey", ls=":", label=f"Brennpunkt F = {F:.1f} mm")
u.set_xlabel("z (mm), Linsenmitte bei 0 (rot: Spalte)")
u.set_ylabel("r (mm)")
u.set_ylim(-1.5, 1.5)
u.set_xlim(-40, max(80, (F or 0) + 20))
u.legend()
fig.tight_layout()
fig.savefig("einzellinse_einfach.png", dpi=130)
print("Bild gespeichert: einzellinse_einfach.png")

# ------------------------------ Kontrollwerte VON CLAUDE via LAPLACE Simulation berechnet. Keine Gewähr! --------------------------------
# Brennpunkt F aus der großen Simulation (Laplace-Gitter, D=10, L=10, g=2,
# 24 eV), zum Vergleich mit eurer Tabelle oben:
#   U = -22: 9,5   -20: 15,0   -18: 22,3   -16: 33,0   -14: 49,1
#   U = -12: 75,3  -10: 120,9   -8: 209,1   +20: 92,6  +30: 53,4
#   U = +45: 33,0  +60: 24,4
# Die Handformel liegt ca. 2-3 % kürzer (bei -22 V ca. 4 %).
# Wenn eure Werte so aussehen, rechnet das Skript richtig.
