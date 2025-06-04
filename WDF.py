import pywdf
from pywdf.core.circuit import Resistor, Capacitor, Inductor, VoltageSource

print(pywdf.__version__)

# === PARAMETRI DI BASE ===
Re = 4.0          # Ω
Cp = 2.4e-8       # Fclear
alpha = 3.7e-4    # N/V
Rm = 9.7e-3       # N·s/m
Mm = 1.0e-6       # kg
Cm = 2.2e-3       # m/N
S_eff = 2.0e-5    # m²
C_bc = 3.6e-3     # Pa/m³
L_tube1 = 1.0e2   # Pa·s²/m³
L_tube2 = 1.0e2   # Pa·s²/m³
C_tube = 6.5e-13  # Pa/m³
R_ac = 5.0e6      # Pa·s/m³

# === CALCOLO VALORI COMPONENTI ===
R1 = Re * alpha**2
R2 = Rm
R3 = R_ac * S_eff**2

L1 = Mm
L2 = L_tube1 * S_eff**2
L3 = L_tube2 * S_eff**2

C1 = Cp / alpha**2
C2 = Cm
C3 = C_bc / S_eff**2
C4 = C_tube / S_eff**2

# === CREAZIONE COMPONENTI WDF ===
Fin = VoltageSource(0)  # sorgente di tensione (input)
r1 = Resistor(R1)
r2 = Resistor(R2)
r3 = Resistor(R3)
l1 = Inductor(L1)
l2 = Inductor(L2)
l3 = Inductor(L3)
c1 = Capacitor(C1)
c2 = Capacitor(C2)
c3 = Capacitor(C3)
c4 = Capacitor(C4)

# === COSTRUZIONE DELLA RETE WDF SECONDO LA TUA DESCRIZIONE ===
a = SeriesAdaptor(Fin, r1)
b = ParallelAdaptor(a, c1)
c = SeriesAdaptor(b, r2)
c = SeriesAdaptor(c, l1)
c = SeriesAdaptor(c, c2)
d = ParallelAdaptor(c, c3)
e = SeriesAdaptor(d, l2)
f = ParallelAdaptor(e, c4)
g = SeriesAdaptor(f, l3)
h = ParallelAdaptor(g, r3)

root = h

# === Visualizzazione struttura ===
print_tree(root)
plot_tree(root)
