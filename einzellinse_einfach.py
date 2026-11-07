r"""
================================================================================
 EINZELLINSE - die einfache Rechnung zum Verstehen
================================================================================

Dieses Skript rechnet aus, wo eine elektrostatische Einzellinse unseren
Elektronenstrahl fokussiert. Es ist absichtlich kurz (der eigentliche Code
sind ~35 Zeilen) und dafür sehr ausführlich kommentiert. Genauigkeit: etwa
2-3 % (bei sehr starker Linse um -22 V ca. 4 %) gegenüber der großen Simulation (Laplace-Gitter + Bahnverfolgung).
Für einen ersten Entwurf reicht das völlig.

Benötigt nur:  pip install numpy matplotlib

--------------------------------------------------------------------------------
 1. POTENTIAL UND ELEKTRONVOLT - was die Zahlen bedeuten
--------------------------------------------------------------------------------
Ein Potential V (in Volt) sagt, wie viel Energie eine Ladung an einem Ort hat,
pro Ladung. Ein Elektron, das von -24 V nach 0 V fliegt, gewinnt
24 Elektronvolt (eV) Bewegungsenergie. Genau das passiert in unserer Kanone:
Kathode auf -24 V, Anode auf 0 V -> die Elektronen kommen mit 24 eV heraus.

Wichtigster Trick des ganzen Skripts: Wir rechnen nicht mit Geschwindigkeiten,
sondern mit der Größe

        Phi(z) = V(z) + 24 V

Phi ist "wie viele Volt Bewegungsenergie das Elektron an der Stelle z gerade
hat". Bei V = 0 (außerhalb der Linse) ist Phi = 24 V, also 24 eV. Liegt die
Mittelelektrode auf -16 V, sinkt Phi dort auf etwa 24 - 16 = 8 V: das
Elektron wird langsamer. Würde Phi irgendwo 0, bliebe das Elektron stehen und
flöge zurück. Deshalb muss |U| kleiner als 24 V bleiben (praktisch <= 0,9*24).

--------------------------------------------------------------------------------
 2. WIE EINE EINZELLINSE FUNKTIONIERT
--------------------------------------------------------------------------------

        Ring 1 (0 V)    Ring 2 (U)     Ring 3 (0 V)
       ===========   =============   ===========      <- Rohrwand (Radius a)
   e- -------->                                    -------> Achse (z)
       ===========   =============   ===========
                  ^                ^
               Spalt 1          Spalt 2

Drei Rohre hintereinander, die äußeren auf 0 V, das mittlere auf U.
In den Spalten "beult" das elektrische Feld ins Rohr hinein. Nahe der Achse
hat dieses Feld eine Komponente, die zur Achse hin oder von ihr weg zeigt,
und zwar umso stärker, je weiter das Elektron von der Achse entfernt ist.
Genau das macht eine Glaslinse mit Licht: weit außen stärker knicken als
innen -> alle Strahlen treffen sich in einem Punkt, dem Brennpunkt.

Jeder Spalt hat eine sammelnde und eine zerstreuende Hälfte. Warum bleibt
insgesamt eine SAMMELNDE Wirkung übrig, egal ob U positiv oder negativ ist?
In der sammelnden Hälfte ist das Elektron langsamer (bzw. weiter außen) und
bleibt länger im Feld als in der zerstreuenden. Unterm Strich gewinnt immer
das Sammeln. "Einzel" heißt die Linse, weil das Elektron hinten mit genau
derselben Energie herauskommt, mit der es hineingeflogen ist.

--------------------------------------------------------------------------------
 3. DIE ACHSNÄHERUNG ("paraxial")
--------------------------------------------------------------------------------
Unser Strahl ist dünn (ca. 1 mm) im Vergleich zur Bohrung (10 mm). Für
Elektronen nahe der Achse muss man das Feld nur AUF der Achse kennen, V(z).
Neben der Achse folgt es automatisch aus der Physik (Laplace-Gleichung):

        V(r, z)  ~  V(z) - r^2/4 * V''(z)

(V'' = zweite Ableitung = "Krümmung" des Potentials entlang z.)
Daraus ergibt sich die Kraft nach außen/innen: proportional zu r. Das
ist der Grund, warum das Ganze überhaupt wie eine Linse wirkt.

--------------------------------------------------------------------------------
 4. DAS POTENTIAL AUF DER ACHSE - die Handformel
--------------------------------------------------------------------------------
Für zwei Rohre (Radius a) mit Spannungssprung dazwischen steigt das Potential
auf der Achse nicht als Treppe an, sondern S-förmig, verschmiert über etwa
einen Rohrradius. Diese S-Kurve ist sehr gut durch tanh beschrieben:

        Sprung(z) = tanh(w * z / a),       w = 1,318  (Erfahrungswert)

Weil unser Spalt eine Breite g hat (die Spannung geht über den Spalt linear
von 0 auf U), mittelt man tanh über die Spaltbreite. Das ergibt:

        S(x) = a/(w*g) * ln( cosh(w*(x+g/2)/a) / cosh(w*(x-g/2)/a) )

S läuft von -1 (weit links) bis +1 (weit rechts). Für die Einzellinse setzt
man einen Sprung hoch (Spalt 1 bei z = -c) und einen runter (Spalt 2 bei +c):

        V(z) = U/2 * ( S(z + c) - S(z - c) ),     c = (L + g)/2

Diese Formel ist nur 1 % (von U) neben der exakten Lösung. Das ist die
Hauptquelle der 2-3 % Fehler im Brennpunkt.

--------------------------------------------------------------------------------
 5. DIE BAHNGLEICHUNG - das Herz der Rechnung
--------------------------------------------------------------------------------
r(z) = Abstand des Elektrons von der Achse an der Stelle z.
r'   = Steigung der Bahn (wie schräg das Elektron fliegt),  r'' = Krümmung.

        r'' = - Phi'/(2*Phi) * r'  -  Phi''/(4*Phi) * r
                 \_____________/       \______________/
                     Term A                Term B

Term B: die eigentliche Linsenkraft. Proportional zu r (s. Abschnitt 3),
        also: weiter außen -> stärker zur Achse gebogen.
Term A: wird das Elektron schneller (Phi' > 0), kommt es in Flugrichtung
        schneller voran als zur Seite -> die Bahn wird "flacher".
        Beim Abbremsen umgekehrt. Ändert also nur die Steigung.

Bemerkenswert: Elektronenmasse und -ladung kommen gar nicht vor! Nur die
Spannungen. Deshalb skaliert alles so schön (Abschnitt 8).

Für Neugierige, woher die Gleichung kommt (Newton + Energieerhaltung):
  (a) Kraft zur Achse:  m * d²r/dt² = -e * r/2 * V''      (aus Abschnitt 3)
  (b) Geschwindigkeit längs der Achse aus der Energie:
          1/2 m v² = e * Phi   ->   v = sqrt(2 e Phi / m)
  (c) Statt "pro Zeit" wollen wir "pro Millimeter z":  dr/dt = r' * v.
      Noch einmal ableiten, (b) einsetzen, durch e teilen -> die Gleichung
      oben. Masse m und Ladung e kürzen sich dabei heraus.

--------------------------------------------------------------------------------
 6. WIE DER COMPUTER DIE GLEICHUNG LÖST - Schritt für Schritt
--------------------------------------------------------------------------------
Wie ein Daumenkino: Wir starten weit vor der Linse mit
r = 1 mm (Abstand von der Achse) und r' = 0 (fliegt parallel zur Achse).
Dann gehen wir in winzigen Schritten dz = 0,01 mm vorwärts und machen
jedes Mal dasselbe:
    1. Krümmung r'' aus der Gleichung ausrechnen (mit Phi an dieser Stelle)
    2. Steigung anpassen:   r' = r' + r'' * dz
    3. Position anpassen:   r  = r  + r'  * dz
Das sind 31.000 Schritte, die der Rechner in Bruchteilen einer Sekunde macht.
(Das Verfahren heißt "semi-implizites Euler-Verfahren". Kleiner dz = genauer.)

--------------------------------------------------------------------------------
 7. WAS AM ENDE HERAUSKOMMT
--------------------------------------------------------------------------------

   r                    Linse
   |  ------------------\                      paralleler Strahl kommt rein
   |                     \__
   |                        \___
   |                            \___
   +------------------|-----|-------X-------> z
                      0     H       F
                            |<--f-->|

  F  Brennpunkt: wo das Elektron die Achse kreuzt (in mm hinter Linsenmitte)
  f  Brennweite: aus der Steigung hinter der Linse, f = -r0 / r'
  H  Hauptebene: von hier aus wird f gemessen, H = F - f. Liegt bei uns fast
     in der Linsenmitte (wenige mm Abweichung), deshalb F ~ f.

Mit f kann man dann wie in der Schul-Optik rechnen (Linsengleichung):
        1/g + 1/b = 1/f      g = Abstand Crossover der Kanone -> H
                             b = Abstand H -> Bild (Fokus)
Das brauchen wir, weil der echte Strahl aus der Kanone nicht parallel ist,
sondern von einem Punkt (Crossover) auseinanderläuft.

--------------------------------------------------------------------------------
 8. SKALIERUNG - die wichtigste Erkenntnis zum Auslegen
--------------------------------------------------------------------------------
* Alle Längen verdoppeln (Bohrung, Ringe, Spalt)  -> Brennpunkt doppelt so weit.
* Strahlenergie UND Linsenspannung verdoppeln     -> Brennpunkt bleibt gleich.
  Es zählt nur das Verhältnis U / E_Strahl.
Ein 48-eV-Strahl mit -32 V verhält sich exakt wie 24 eV mit -16 V.

--------------------------------------------------------------------------------
 9. WAS DIESE RECHNUNG NICHT KANN (bewusst weggelassen)
--------------------------------------------------------------------------------
* Elektronen weit außen (> ca. 1/3 des Radius): dort wird der Fokus kürzer
  (Öffnungsfehler). Das kann nur die große Simulation.
* Energieunschärfe: 24 +- 0,5 eV verschiebt den Fokus bei -16 V um +-2 mm.
* Raumladung (die Elektronen stoßen sich gegenseitig ab): bei 24 eV schon ab
  ~1 µA merklich -> Strahlstrom im nA-Bereich halten.
* Erdmagnetfeld: lenkt 24-eV-Elektronen auf 33 mm schon ~1,5 mm ab.
================================================================================
"""
import numpy as np
import matplotlib.pyplot as plt

# =========================== EINGABEN (mm, V, eV) ============================
D = 10.0        # Bohrung = Innendurchmesser der Ringe
L = 10.0        # Länge des MITTLEREN Rings (die äußeren sind hier egal,
                # solange sie mindestens ~D lang sind)
g = 2.0         # Spaltbreite zwischen den Ringen
E = 24.0        # Strahlenergie in eV (= Kathodenspannung, nur ohne Minus)
U_BILD = -16.0  # für diese Spannung wird ein Bild gezeichnet
SPANNUNGEN = [-22, -20, -18, -16, -14, -12, -10, -8, 20, 30, 45, 60]
# =============================================================================

a = D / 2                       # Rohrradius
w = 1.318                       # Formfaktor der tanh-Näherung (Abschnitt 4)
c = (L + g) / 2                 # Abstand Linsenmitte -> Spaltmitte
dz = 0.01                       # Schrittweite in mm (Abschnitt 6)
z = np.arange(-60, 250, dz)     # alle Stellen z, an denen wir rechnen.
                                # Linsenmitte bei z = 0, Start 60 mm davor.


def S(x):
    """Ein S-förmiger Sprung von -1 auf +1 über einen Spalt (Abschnitt 4)."""
    return a / (w * g) * np.log(np.cosh(w * (x + g/2) / a) / np.cosh(w * (x - g/2) / a))


def bahn(U, r0=1.0):
    """Verfolgt ein Elektron, das im Abstand r0 parallel zur Achse startet.
    Gibt das Potential V(z), die Bahn r(z), den Brennpunkt F, die
    Brennweite f und die Hauptebene H zurück."""
    V = U / 2 * (S(z + c) - S(z - c))      # Potential auf der Achse
    Phi = V + E                            # Bewegungsenergie in Volt
    if Phi.min() <= 0:                     # Elektron bleibt stehen -> Spiegel
        return V, None, None, None, None
    Phi1 = np.gradient(Phi, dz)            # Phi'  (Steigung von Phi)
    Phi2 = np.gradient(Phi1, dz)           # Phi'' (Krümmung von Phi)

    r, rs = r0, 0.0                        # Startabstand, Startsteigung r'
    rr = np.empty_like(z)                  # hier merken wir uns die Bahn
    for i in range(len(z)):                # das Daumenkino (Abschnitt 6)
        rr[i] = r
        r2 = -Phi1[i] / (2 * Phi[i]) * rs - Phi2[i] / (4 * Phi[i]) * r
        rs += r2 * dz                      # Steigung anpassen
        r += rs * dz                       # Position anpassen
        if z[i] >= 60 and i > 0 and z[i - 1] < 60:
            r60, rs60 = r, rs              # Werte weit hinter der Linse merken

    # Brennpunkt F: erste Stelle hinter der Mitte, an der r von + nach - geht
    kreuz = np.where((z[1:] > 0) & (rr[:-1] > 0) & (rr[1:] <= 0))[0]
    F = None
    if len(kreuz):
        k = kreuz[0]                       # zwischen zwei Punkten linear
        F = z[k] + rr[k] / (rr[k] - rr[k + 1]) * dz   # interpolieren
    f = -r0 / rs60                         # Brennweite aus der Endsteigung
    H = (60 - r60 / rs60) - f              # Hauptebene = F(Asymptote) - f
    return V, rr, F, f, H


# ------------------------- Tabelle für alle Spannungen -----------------------
print(f"Einzellinse: Bohrung {D} mm, Mittelring {L} mm, Spalt {g} mm, Strahl {E} eV")
print("  U (V) |  U/E  | Brennpunkt F | Brennweite f | Hauptebene H   (mm, ab Linsenmitte)")
for U in SPANNUNGEN:
    V, rr, F, f, H = bahn(U)
    if rr is None:
        print(f"{U:7.0f} | {U/E:+.2f} | Elektron wird reflektiert (|U| >= E)")
    else:
        Ftxt = f"{F:8.1f}" if F else "  > 250 "
        print(f"{U:7.0f} | {U/E:+.2f} |   {Ftxt}   |   {f:8.1f}   |   {H:+6.1f}")

# ------------------------------- Bild zeichnen -------------------------------
V, rr, F, f, H = bahn(U_BILD)
fig, (o, u) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
o.plot(z, V + E, color="#4a3aa7")
o.axhline(E, color="grey", lw=0.5, ls="--")
o.set_ylabel("Phi = Bewegungsenergie (eV)")
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

# ------------------------------ Kontrollwerte --------------------------------
# Brennpunkt F aus der großen Simulation (Laplace-Gitter, D=10, L=10, g=2,
# 24 eV), zum Vergleich mit eurer Tabelle oben:
#   U = -22: 9,5   -20: 15,0   -18: 22,3   -16: 33,0   -14: 49,1
#   U = -12: 75,3  -10: 120,9   -8: 209,1   +20: 92,6  +30: 53,4
#   U = +45: 33,0  +60: 24,4
# Die Handformel liegt ca. 2-3 % kürzer (bei -22 V ca. 4 %). Wenn eure Werte so aussehen,
# rechnet das Skript richtig.
