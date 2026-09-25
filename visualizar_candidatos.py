import matplotlib.pyplot as plt
import numpy as np

# Datos de las 2 variantes candidatas en BUB1B (de candidatos_mva.py + VEP)
# Color = consecuencia funcional (rojo = stop_gained/HIGH, naranja = missense/MODERATE)
# Gris = referencia / poblacion general, para que el ojo distinga "hallazgo" de "fondo"
COLOR_STOP = "#e34948"      # rojo (stop_gained, HIGH impact)
COLOR_MISSENSE = "#2a78d6"  # azul (missense, MODERATE impact)
COLOR_NEUTRO = "#c3c2b7"    # gris neutro (referencia / fondo poblacional)

variantes = [
    {
        "etiqueta": "chr15:40209701 T>G\np.Leu737Ter (stop_gained)",
        "color": COLOR_STOP,
        "pos_aa": 737,
        "ad_ref": 21,
        "ad_alt": 25,
        "gnomad_af": 3.286e-05,
        "af_observada": True,
    },
    {
        "etiqueta": "chr15:40220612 T>G\np.Asn1002Lys (missense)",
        "color": COLOR_MISSENSE,
        "pos_aa": 1002,
        "ad_ref": 15,
        "ad_alt": 13,
        "gnomad_af": 1e-6,  # piso visual: no observada en gnomAD
        "af_observada": False,
    },
]
LONGITUD_PROTEINA_BUB1B = 1050  # aa, ENSP00000287598

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# --- Panel A: profundidad alelica (evidencia de que el genotipo 0/1 es solido) ---
y = np.arange(len(variantes))
alto = 0.32
for i, v in enumerate(variantes):
    axes[0].barh(i + alto / 2, v["ad_ref"], height=alto, color=COLOR_NEUTRO, label="Reads REF" if i == 0 else None)
    axes[0].barh(i - alto / 2, v["ad_alt"], height=alto, color=v["color"], label="Reads ALT" if i == 0 else None)
    axes[0].text(v["ad_ref"] + 0.5, i + alto / 2, str(v["ad_ref"]), va="center", fontsize=9)
    axes[0].text(v["ad_alt"] + 0.5, i - alto / 2, str(v["ad_alt"]), va="center", fontsize=9)

axes[0].set_yticks(y)
axes[0].set_yticklabels([v["etiqueta"] for v in variantes], fontsize=9)
axes[0].set_xlabel("Numero de reads")
axes[0].set_title("Profundidad alelica (AD)\nbalance ~50/50 = genotipo 0/1 confiable")
axes[0].legend(loc="lower right", frameon=False)
axes[0].spines[["top", "right"]].set_visible(False)

# --- Panel B: rareza poblacional (paciente vs gnomAD, escala log) ---
y2 = np.arange(len(variantes))
for i, v in enumerate(variantes):
    axes[1].barh(i + alto / 2, 0.5, height=alto, color=v["color"], label="Paciente (heterocigota)" if i == 0 else None)
    axes[1].barh(i - alto / 2, v["gnomad_af"], height=alto, color=COLOR_NEUTRO,
                 hatch="" if v["af_observada"] else "//",
                 label="gnomAD (poblacion general)" if i == 0 else None)
    texto_gnomad = f"{v['gnomad_af']:.1e}" if v["af_observada"] else "no observada"
    axes[1].text(max(v["gnomad_af"], 1e-6) * 1.3, i - alto / 2, texto_gnomad, va="center", fontsize=9)
    axes[1].text(0.5 * 1.15, i + alto / 2, "0.5", va="center", fontsize=9)

axes[1].set_xscale("log")
axes[1].set_xlim(1e-6, 1)
axes[1].set_yticks(y2)
axes[1].set_yticklabels([v["etiqueta"] for v in variantes], fontsize=9)
axes[1].set_xlabel("Frecuencia alelica (escala log)")
axes[1].set_title("Rareza poblacional\npaciente vs. gnomAD")
axes[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=False, fontsize=9)
axes[1].spines[["top", "right"]].set_visible(False)

# --- Panel C: ubicacion en la proteina BUB1B (needle/lollipop plot) ---
axes[2].hlines(0, 0, LONGITUD_PROTEINA_BUB1B, color=COLOR_NEUTRO, linewidth=6, zorder=1)
for v in variantes:
    axes[2].vlines(v["pos_aa"], 0, 1, color=v["color"], linewidth=2, zorder=2)
    axes[2].plot(v["pos_aa"], 1, "o", color=v["color"], markersize=10, zorder=3)
    axes[2].annotate(
        v["etiqueta"].split("\n")[1],
        xy=(v["pos_aa"], 1),
        xytext=(v["pos_aa"], 1.35),
        ha="center",
        fontsize=8,
        rotation=15,
    )

axes[2].set_xlim(0, LONGITUD_PROTEINA_BUB1B)
axes[2].set_ylim(-0.3, 2.1)
axes[2].set_yticks([])
axes[2].set_xlabel("Posicion en la proteina BUB1B (aa)")
axes[2].set_title("Localizacion de las variantes\nsobre la proteina (1050 aa)")
axes[2].spines[["top", "right", "left"]].set_visible(False)

fig.suptitle("Variantes candidatas en BUB1B (heterocigoto compuesto propuesto)", fontsize=13, y=1.05)
plt.tight_layout()
fig.subplots_adjust(top=0.78, bottom=0.22)
plt.savefig("variantes_candidatas_bub1b.png", dpi=150, bbox_inches="tight")
plt.show()
