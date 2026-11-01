# NOVELTY_AUDIT — Estado del arte y novedad de `DISCOVERY_REPORT.md`
*Generado por `09_make_audit.py` desde `discovery_results.json` (incl. `08_novelty_checks.py`). 2026-11-01.*

> **Actualización tras los cierres teóricos:** P10/P11 (Penrose cargada derivada internamente, con alcance acotado) y 05b (cinética con bloqueo de Pauli). L1 sigue siendo Cat.1; K1 se añade como Cat.2 dependiente de modelo.

## 0. Alcance y límites (léase primero)
- **Ningún texto completo fue leído.** `arxiv.org` y `ar5iv.labs.arxiv.org` no resuelven desde esta sesión (red restringida): hay que permitir esos dominios en *Network access* del entorno para una auditoría de texto completo. Todo lo que sigue se apoya en **extractos de búsqueda** (nivel **E**) o en **memoria** (nivel **M**, no verificada). La *ausencia* de un resultado en una búsqueda superficial (nivel **N**) **no** prueba que no exista.
- Esta auditoría **no** puede dar «rigor absoluto» sobre la novedad: da una clasificación conservadora y una lista de verificaciones pendientes (§7).
- **Corrección de resultados propios** hallada durante la auditoría (ya aplicada en 04, informe y JSON): el margen DEC «−8.33e-4, acotado» era el borde del rango (ε=10⁴⁰); el umbral «1.92» era artefacto de rejilla (exacto **1.88218**). Ahora hay 67 gates.

## 1. Veredicto en una tabla
| Resultado | Categoría | En una línea |
|---|---|---|
| **L1** Cancelación de Coulomb en el atrapamiento de un cascarón (P1–P3) | Cat.1 Redescubrimiento | G·(M−Q²/8πR)=R/2 ⇔ g^{rr}(R)=0 (definición del horizonte de RN). |
| **L2a** Escala τ_dis ∝ (eE₀)^{-1/2} (exponente −0.5000) | Cat.1 Redescubrimiento | t·√(eE₀)=κ(e) por invariancia de escala. |
| **L2b** Constante κ = 13.5394 | Cat.2 Refinamiento | Valor numérico con convergencia (orden 1.98) y método independiente. |
| **L3** Ley τ_dis/τ_ff=√(r₀*/r₀), r₀*=9.29e-18 l_P, independiente de N | Cat.3 Síntesis | τ_dis/τ_ff ∼ (Gm²/e²)^{1/2}·(κ-factores); z²=4.17e+42. |
| **L4a** Crecimiento logarítmico a un bucle y polo de Landau | Cat.1 Redescubrimiento | L=(E²/2)[1−(α/3π)(ln ε+const)]. |
| **L4b** Violación de DEC desde ε=1.88218; margen -7.743e-04 (asintótico) | Cat.1 Redescubrimiento | margen→−α/(3π)=-7.743e-04; a ε=10¹⁰: -7.875e-04. |
| **A1** Perfil interior estático: EH cuártico vs un bucle exacto vs Born–Infeld (04) | Cat.2 Refinamiento | E∼r^{-2/3} y masa finita del cuártico como artefacto (S9); un bucle exacto = RN con carga corrida (δr₊/r₊ medido -1.191e-02 vs predicho -1.230e-02); BI/cuártico → horizonte único tipo Schwarzschild. |
| **B1** Modelo cinético 05 (fuente, conducción + polarización, balance de energía) | Cat.1 Redescubrimiento | Mismas ecuaciones esenciales. |
| **K1** Cinética de Markov con bloqueo de Pauli a través de ≥5 oscilaciones (05b): ω_p=0.1865 e√ε₀, γ/ω_p=0.0052 | Cat.2 Refinamiento | Ecuación de Vlasov cuántica en el límite de Markov con bloqueo de Pauli. |
| **P1** Diagrama de fases (N, r₀) | Cat.3 Síntesis | Colocarlas en un solo plano con valores explícitos, tiempos y verificación simbólica de cada umbral (P1–P9). |
| **T1** Umbrales N*=2.74e+05, N_crit=5.71e+44, N_valid=9/α | Cat.1 Redescubrimiento | N* como cantidad de modelo (Cat.3 menor). |
| **M1** Metodología: pipeline «LLM propone, CAS ejecuta», gates simbólicos, controles negativos, informe autogenerado | Cat.3 Síntesis | 67 gates (15 controles negativos), trazabilidad número→clave JSON. |

Recuento: Cat.1 = 6, Cat.2 = 3, Cat.3 = 3. **Ningún resultado es, hoy, un descubrimiento físico nuevo defendible**; el valor está en la síntesis cuantitativa, la corrección explícita de la truncación y la verificación formal.

## 2. Comparación artículo por artículo (nivel de verificación entre paréntesis)
| Trabajo | Nivel | Qué se sabe (según lo visto) | Gate(s) del repositorio | Relación |
|---|---|---|---|---|
| Heisenberg & Euler, Z. Phys. 98, 714 (1936); Weisskopf, Kgl. Danske Vid. Selsk. Mat.-fys. Medd. 14 (6) (1936) | M | Lagrangiano efectivo de QED a un bucle; término cuártico (2α²/45m⁴)[(E²−B²)²+7(E·B)²]. | 04:S5,S6 (coef. 4cE⁴), 01/03 (baseline cuártico) | Cat.1 Redescubrimiento |
| Schwinger, Phys. Rev. 82, 664 (1951) | M | Forma de tiempo propio exacta; Im L = tasa de decaimiento del vacío; serie en e^{−nπ/ε}. | 04:N1 (Im L₁ = serie, 8e-15), S6; 05:D3 (tasa de pares) | Cat.1 Redescubrimiento |
| Ritus, Zh. Eksp. Teor. Fiz. 69, 1517 (1975); Ritus hep-th/9812124 (1998); Lebedev & Ritus (representación dispersiva) | E | Aparece como referencia estándar del Lagrangiano de campo intenso y QED a distancias cortas; un extracto de 2024 indica que el límite de campo fuerte es logarítmico con el coeficiente de la función β (Lebedev–Ritus: representación dispersiva). No se vio el texto. | 04:S7 (coef. −α/3π), S8 (polo), N3 (pendiente −1/24π²), N8 | Cat.1 Redescubrimiento |
| Gibbons, Commun. Math. Phys. 44, 245 (1975) | E | Pérdida espontánea de carga de agujeros negros por emisión de partículas; para agujeros grandes (>~10¹⁷ g) rige una fórmula tipo Schwinger. Un trabajo de 2025 (arXiv:2512.01898) resume que Gibbons calculó la probabilidad de descarga con coeficientes de Bogoliubov para un campo escalar cargado, que los agujeros pequeños (campo local mayor) se descargan más, y que esto implica una cota inferior de masa para agujeros extremales. | 06:P5 (N_crit=M_Pl²/m²), 04 (ε(r_Q)=M_Pl²/N), §7 «facilitador» | Cat.1 Redescubrimiento |
| Drummond & Hathrell, Phys. Rev. D 22, 343 (1980) | E | Según trabajos posteriores: acción efectiva a un bucle con términos RF², R_μνF², R_μνλρF² y un término derivativo, con coeficientes O(α)/m²; desplazamiento de velocidad de fotones O(α(λ_C/L)²), nulo en horizontes. | NINGUNO: nuestro tratamiento curvo es de acoplo mínimo. LAGUNA: a L≈r_Q≪λ_C el cociente λ_C²/L² es enorme y la expansión no está controlada | Laguna / tensión |
| Yajima & Tamaki, Phys. Rev. D 63, 064007 (2001), gr-qc/0005016 | E | Agujeros negros estáticos esféricos en Einstein–Euler–Heisenberg (EH como teoría efectiva), caso magnético resuelto analíticamente; parámetros a,b; dependencia de causalidad y termodinámica; singularidad central; a>0 natural si EH es límite de baja energía de Born–Infeld. | 01/03 (EEH cuártico estático), 04:S5 (relación BI→EH) | Cat.1 Redescubrimiento |
| Ruffini, Wu & Xue, arXiv:1307.4951 (2013) [revista: Phys. Rev. D 88, 085004 – M] | E | Einstein + Lagrangiano EH a un bucle; casos eléctrico, magnético y diónico; horizonte, área, entropía, energía total y extraíble con corrección QED; la condición de agujero extremo se desplaza. Resultados a orden dominante en la corrección QED; el desarrollo O(α³) se escribe. Una búsqueda de seguimiento no encontró cálculo no perturbativo (tiempo propio exacto) de horizonte RN. | 01 (derivación), 03/04 (cuártico); 04 extiende a tiempo propio exacto (Cat.2) | Cat.1 Redescubrimiento |
| Cooper & Mottola, Phys. Rev. D 40, 456 (1989) [volumen/página: M] | M | Retroacción cuántica como problema de valores iniciales (QED escalar). No localizado en las búsquedas. | 05 (contexto) | Cat.1 Redescubrimiento |
| Kluger, Eisenberg, Svetitsky, Cooper, Mottola, PRL 67, 2427 (1991); PRD 45, 4659 (1992) [M]; Int. J. Mod. Phys. E (1993) = hep-ph/0311293 | E | Producción de pares en campo eléctrico fuerte como problema de valores iniciales en campo medio semiclásico homogéneo; la retroacción genera una corriente que se opone al campo y aparecen oscilaciones de plasma. Se compara con un modelo Boltzmann–Vlasov con fuente de Schwinger y factor de bloqueo de Pauli (1−2f). | 05 (modelo completo: fuente, corriente de conducción + polarización, D2 balance de energía, oscilaciones) | Cat.1 Redescubrimiento |
| Kluger, Mottola & Eisenberg, hep-ph/9803372 (1998) [PRD 58, 125015 – M] | E | El número de partículas adiabático obedece una ecuación de Vlasov cuántica no local en el tiempo; para campos débiles y lentos actúa como función de distribución y su derivada como fuente de Boltzmann–Vlasov (límite de Markov). | 05 (fuente local = límite de Markov; nosotros: solo término n=1, sin (1−2f), válido hasta el primer cero) | Cat.1 Redescubrimiento |
| Schmidt et al., hep-ph/9809227, Int. J. Mod. Phys. E 7 (1998) | E | Ecuación cinética cuántica no markoviana; a baja densidad el término fuente se reduce a la parte principal de la fórmula de Schwinger. | 05 (justificación de la fuente local) | Cat.1 Redescubrimiento |
| Ruffini, Vereshchagin & Xue, Phys. Rep. 487, 1 (2010), arXiv:0910.0974 | E | Revisión; §8 trata oscilaciones de plasma, retroacción y apantallamiento por pares. Un extracto (de segunda mano) da longitudes de oscilación: E₀=10E_c → ~10 ħ/m_ec; E₀=0.15E_c → ~10⁷ ħ/m_ec. | 05: contraste cuantitativo (08): nuestro primer cero 17.1 y 8.78e+06 (cociente 5.1e+05 vs 10⁶ citado) | Cat.1 Redescubrimiento |
| Gelis & Tanji, Prog. Part. Nucl. Phys. 87, 1 (2016), arXiv:1510.05451 | E | Solo metadatos y resumen (revisión del mecanismo de Schwinger, incl. Schwinger asistido dinámicamente). Contenido NO accedido. | 05 (no comparado; pendiente) | Laguna / tensión |
| Prakapenia & Vereshchagin, arXiv:2306.12956 (2023) | E | El bloqueo de Pauli suprime la creación de pares en campos sobrecríticos y modifica la retroacción; cuando los pares dominan energéticamente la frecuencia de oscilación aumenta. | 05: límite declarado (validez hasta el primer cero) | Cat.1 Redescubrimiento |
| Boulware, Phys. Rev. D 8, 2363 (1973) | E | Fuente secundaria: capa fina cargada con uniones de Israel y extensión máxima de RN; ecuación de energía que relaciona masa total, masa propia, carga y radio; singularidad desnuda solo con masa propia negativa. | 06:P1–P3 (L1), §7 | Cat.1 Redescubrimiento |
| Bekenstein, Phys. Rev. D 4, 2185 (1971) | E | Equilibrio hidrostático y colapso de bolas de fluido cargado; se cita como origen del «doble papel» de la carga y del teorema de Birkhoff cargado. Contenido detallado NO accedido. | 06 (L1, z²), §7 | Cat.1 Redescubrimiento |
| Malec & Ó Murchadha, Phys. Rev. D 49 (1994) [M], gr-qc/9401014 | E | Desigualdad de Penrose en simetría esférica aun con el horizonte dentro de la materia; condición necesaria para formación de superficies atrapadas; la versión cargada tipo Gibbons se demuestra si toda la carga está dentro de la superficie atrapada exterior y, según el resumen, puede romperse en etapas tempranas del colapso. | 06:P1–P3, criterio (ii); riesgo: hipótesis de carga contenida | Cat.1 Redescubrimiento |
| Khuri, Weinstein & Yamada: arXiv:1308.3771 (2013), arXiv:1409.3271 (2014/15), Contemp. Math. (2015) | E | Desigualdad de Penrose riemanniana con carga (cota superior r ≤ m+√(m²−q²)) bajo DEC cargada y SIN materia cargada fuera del horizonte; m≥|q|. El «2017» del encargo corresponde, de memoria (NO verificado), a la versión publicada de arXiv:1409.3271 en J. Differential Geom. 106 (2017); lo visto en extractos es 2013–2015 (preprints/Contemp. Math.). | 06:P1 (RN satura), criterio (ii) | Cat.1 Redescubrimiento |
| Misner & Sharp, Phys. Rev. 136, B571 (1964) [M]; Hayward (1994) trapping horizons [E parcial]; Hod (masa M∞−Q²/2R) [E vía arXiv:2606.31116] | E | Un trabajo reciente indica que, en RN, la masa dentro de R, M(r≤R)=M∞−Q²/2R, coincide exactamente con la masa cuasilocal de Misner–Sharp M=(R/2)(1−g^{rr}); el horizonte de atrapamiento es el borde 2M(R)=R. | 06:P2, P3 = ESTA identidad. L1 es una reescritura. | Cat.1 Redescubrimiento |
| Russo & Townsend, JHEP 06 (2024) 191, arXiv:2404.09994 | E | En NLED de clase Plebanski, causalidad implica DEC y SEC; toda violación implica acausalidad. Tomizawa & Suzuki (arXiv:2309.10535): violación de DEC ↔ fotones superlumínicos en la geometría efectiva. | 04:N8, N9, S10 (nuestra DEC solo en configuración eléctrica pura, P=0): TENSIÓN a comprobar | Laguna / tensión |
| Fedotov et al., arXiv:1608.02261; Mironov, Meuren & Fedotov (2020) [identificador arXiv:2003.06909 NO confirmado]; Heinzl, Ilderton & King (unitariedad) | E | Conjetura de Ritus–Narozhny: cada lazo aporta αχ^{2/3}, no α; la teoría perturbativa fallaría para αχ^{2/3}≳1 (χ≳α^(-3/2)≈1604) en campo cruzado constante. La resumación indica ruptura ya en g≳1; Heinzl–Ilderton–King discuten que la unitariedad elimina diagramas supuestamente responsables. Abierta; se enuncia para campo cruzado constante (aplicabilidad a E puro: NO verificada). | 04:N7,N8; 05 barrido hasta ε=10¹⁰; DISCOVERY_REPORT §4/§6. LAGUNA importante | Laguna / tensión |
| Gershtein & Zel'dovich (1969); Pieper & Greiner (1969); Greiner, Müller & Rafelski (1985) [M]; Popov et al. arXiv:2008.05005, 1910.01373 [E] | E | Núcleo extendido: el nivel 1s alcanza el continuo negativo en Z≈173; carga puntual: el nivel existe solo hasta Z≈137; decaimiento del vacío con emisión espontánea de positrones; apantallamiento y bloqueo de Pauli reducen la producción. | 06: N_valid=9/α, «r_neutral», supercriticidad Nα≳1 | Cat.1 Redescubrimiento |
| Adler, Collins & Duncan, Phys. Rev. D 15, 1712 (1977); Collins, Duncan & Joglekar, Phys. Rev. D 16, 438 (1977) [M] | M | Anomalía de traza T^μ_μ=(β/2e)F². En nuestro resultado: margen DEC asintótico −α/(3π)=-7.743e-04. | 04:S10, N9 | Cat.1 Redescubrimiento |
| Worldline-instanton / Schwinger en RN extremal: arXiv:2409.02197; revisita de Gibbons: arXiv:2512.01898; arXiv:1712.01302 | E | La tasa de Schwinger en RN extremal se apaga al acercarse la razón q/m a la cota de extremalidad; el efecto es máximo en el horizonte y decae en una escala Compton; un RN se descarga apreciablemente cuando m² ≲ (carga·campo)/(radio transversal)². | 06 (N_crit), §7 | Cat.1 Redescubrimiento |

## 3. Clasificación detallada
### L1 — Cancelación de Coulomb en el atrapamiento de un cascarón (P1–P3)  →  **Cat.1 Redescubrimiento**
- **Estructura ya publicada:** Misner–Sharp en RN: M(R)=M∞−Q²/2R (Hod; Misner & Sharp 1964); ecuación de energía de capa fina (Israel; Boulware 1973); desigualdad de Penrose cargada saturada por RN (Gibbons 1972 [M]; Malec–Ó Murchadha 1994; Khuri–Weinstein–Yamada).
- **Ecuación equivalente:** G·(M−Q²/8πR)=R/2 ⇔ g^{rr}(R)=0 (definición del horizonte de RN).
- **Aporte del repositorio:** CIERRE (P10, P11): derivación interna con SymPy de θ₊=0 ⇔ 2GM_MS=R, M_ADM=M_MS+U_C y la igualdad de Penrose cargada para datos simétricos en el tiempo, sin teoremas externos; U_C exacta a un bucle (δU/U_C=7.047e-02 en R=r_Q, N=1e6). Sigue siendo Cat.1: el contenido es conocido; lo que cambia es que ya no depende de una cita.
- **No presentar como:** «Descubrimos que la energía electrostática no cuenta para atrapar».
- **Cómo citar:** «Recordamos que, en términos de la masa de Misner–Sharp M(R)=M−Q²/(2R), el criterio de atrapamiento es 2GM(R)≥R, de modo que la energía electrostática exterior no contribuye [Misner–Sharp 1964; Hod; Boulware 1973].»
- **Salvedades:** Probado internamente sólo para datos esféricos simétricos en el tiempo, ρ≥0 y sin materia cargada exterior (las mismas hipótesis de KWY); el caso dinámico general (Misner–Sharp completo) NO. Malec–Ó Murchadha indican que la versión cargada puede romperse en etapas tempranas del colapso si hay carga exterior.

### L2a — Escala τ_dis ∝ (eE₀)^{-1/2} (exponente −0.5000)  →  **Cat.1 Redescubrimiento**
- **Estructura ya publicada:** Análisis dimensional con m→0 (única escala √(eE)); frecuencia de plasma ω_p∝e√(eE) en la retroacción de Schwinger (Kluger et al. 1991–1998; Ruffini–Vereshchagin–Xue 2010 §8) [fórmula exacta: verificar en texto completo].
- **Ecuación equivalente:** t·√(eE₀)=κ(e) por invariancia de escala.
- **Aporte del repositorio:** Medición numérica del exponente y verificación de convergencia.
- **No presentar como:** «Ley nueva τ∝E^{-1/2}».
- **Cómo citar:** «Consistente con la escala de plasma conocida, ω_p∼e√(eE) [Kluger et al.].»
- **Salvedades:** Exponente = consecuencia dimensional, no resultado.

### L2b — Constante κ = 13.5394  →  **Cat.2 Refinamiento**
- **Estructura ya publicada:** Constante numérica de un modelo concreto (fuente local n=1, LCFA, primer cero del campo, 8 nodos Gauss–Laguerre en p⊥). Contraste con segunda mano: primeros ceros 17.1 (E₀=10E_c) y 8.78e+06 (0.15E_c) frente a ≈10 y ≈10⁷ en un extracto de Ruffini et al.; coinciden en orden de magnitud en ambos extremos (cociente 5.1e+05 vs 10⁶).
- **Ecuación equivalente:** —
- **Aporte del repositorio:** Valor numérico con convergencia (orden 1.98) y método independiente.
- **No presentar como:** «κ=13.5394 es una constante universal de la descarga de Schwinger».
- **Cómo citar:** «En un modelo de Markov con LCFA y primer cero del campo como definición, obtenemos κ≈13.54.»
- **Salvedades:** Dependencia de modelo O(1); sin (1−2f) ni QKE; definición «primer cero» no estándar; comparación incompleta con la literatura.

### L3 — Ley τ_dis/τ_ff=√(r₀*/r₀), r₀*=9.29e-18 l_P, independiente de N  →  **Cat.3 Síntesis**
- **Estructura ya publicada:** Es el producto de dos leyes de potencia con el mismo √N (τ_dis∝r₀/√N, τ_ff∝√(r₀³/N)); r₀* está fijada por e²/(Gm²)=αM_Pl² (el cociente EM/gravedad del electrón, conocido).
- **Ecuación equivalente:** τ_dis/τ_ff ∼ (Gm²/e²)^{1/2}·(κ-factores); z²=4.17e+42.
- **Aporte del repositorio:** Prefactor cerrado con κ y verificación simbólica (P7–P9).
- **No presentar como:** «Ley universal sin precedentes»; r₀*≪l_P está fuera de la validez de ambos modelos (solo la desigualdad τ_dis<τ_ff para r₀>l_P tiene sentido).
- **Cómo citar:** «Para N electrones, el cociente entre el tiempo de descarga y el de caída libre escala como (r₀*/r₀)^{1/2} con r₀* determinado por αM_Pl² y κ.»
- **Salvedades:** A r₀≲λ_C la caída libre de polvo de N masas en reposo no describe materia degenerada confinada.

### L4a — Crecimiento logarítmico a un bucle y polo de Landau  →  **Cat.1 Redescubrimiento**
- **Estructura ya publicada:** Weisskopf 1936; Schwinger 1951; Ritus 1975; Lebedev–Ritus; coeficiente de la función β. Cross-check cuantitativo: nuestro ε_L → μ_L≈10^277.7 GeV frente a ≈10²⁷⁷ GeV citado para el polo de Landau de la QED (extracto de arXiv:1812.08673).
- **Ecuación equivalente:** L=(E²/2)[1−(α/3π)(ln ε+const)].
- **Aporte del repositorio:** Constante asintótica K y evaluación numérica exacta del rayo (Re y Im).
- **No presentar como:** «El crecimiento logarítmico y el polo en log₁₀ε≈562 son resultados nuevos».
- **Cómo citar:** «Reproducimos el límite logarítmico conocido del Lagrangiano de Heisenberg–Euler [Ritus; Lebedev–Ritus; Dunne].»
- **Salvedades:** Solo un bucle. Conjetura de Ritus–Narozhny: αχ^{2/3}≳1 ya en χ≈1604; nuestro barrido llega a ε=10¹⁰ (6.8 órdenes por encima).

### L4b — Violación de DEC desde ε=1.88218; margen -7.743e-04 (asintótico)  →  **Cat.1 Redescubrimiento**
- **Estructura ya publicada:** Anomalía de traza a un bucle (Adler–Collins–Duncan 1977; Collins–Duncan–Joglekar 1977 [M]). Identidad S10: E·D−2L=x³d(L₁/x²)/dx.
- **Ecuación equivalente:** margen→−α/(3π)=-7.743e-04; a ε=10¹⁰: -7.875e-04.
- **Aporte del repositorio:** Umbral exacto ε=1.88218 (máximo de L₁/x²) con integral exacta (Cat.2, menor). CORRECCIÓN: el «−8.33e-4 acotado» del informe previo era el borde del rango (ε=10⁴⁰: -8.328e-04); el margen NO está acotado.
- **No presentar como:** «Violación cuántica acotada de la DEC» como hallazgo; el valor −8.33e-4 como cota.
- **Cómo citar:** «La parte real del Lagrangiano a un bucle viola la DEC en configuración eléctrica pura para ε>1.88; el margen relativo tiende a −α/3π, reflejo de la anomalía de traza.»
- **Salvedades:** Tensión con Russo–Townsend (causalidad⇒DEC): falta calcular las condiciones de causalidad con L₁(S,P) completo; para ε>1 Im L≠0, la parte real aislada no es una acción efectiva estable.

### A1 — Perfil interior estático: EH cuártico vs un bucle exacto vs Born–Infeld (04)  →  **Cat.2 Refinamiento**
- **Estructura ya publicada:** Yajima–Tamaki 2001 y Ruffini–Wu–Xue 2013 trabajan con EH truncado (cuártico/leading order); una búsqueda de seguimiento no halló horizonte de RN con el tiempo propio exacto en campo fuerte (nivel N).
- **Ecuación equivalente:** —
- **Aporte del repositorio:** E∼r^{-2/3} y masa finita del cuártico como artefacto (S9); un bucle exacto = RN con carga corrida (δr₊/r₊ medido -1.191e-02 vs predicho -1.230e-02); BI/cuártico → horizonte único tipo Schwarzschild.
- **No presentar como:** «La QED no censura el horizonte» como predicción física en r_Q.
- **Cómo citar:** «Mostramos que la extrapolación cuártica del EEH produce un perfil cualitativamente distinto del de la integral de tiempo propio exacta; comparamos ambos en el sistema acoplado.»
- **Salvedades:** En los perfiles de 04, la escala de variación del campo es ℓ≈r_Q/2=1.8e-18 λ_C ≪ λ_C: viola el criterio derivativo del propio baseline (03) y la expansión de Drummond–Hathrell no está controlada. Son soluciones FORMALES del modelo de campo localmente constante.

### B1 — Modelo cinético 05 (fuente, conducción + polarización, balance de energía)  →  **Cat.1 Redescubrimiento**
- **Estructura ya publicada:** Kluger et al. 1991–98 (modelo y límite de Markov), Schmidt et al. 1998.
- **Ecuación equivalente:** Mismas ecuaciones esenciales.
- **Aporte del repositorio:** Verificación simbólica D1/D2, implementación Volterra con convergencia de orden 2 y controles negativos.
- **No presentar como:** «Nuevo modelo de retroacción».
- **Cómo citar:** «Reimplementamos el modelo de Markov de Kluger–Mottola–Eisenberg con verificación simbólica del balance de energía.»
- **Salvedades:** Sin (1−2f), sin QKE, sin transporte espacial; válido hasta el primer cero.

### K1 — Cinética de Markov con bloqueo de Pauli a través de ≥5 oscilaciones (05b): ω_p=0.1865 e√ε₀, γ/ω_p=0.0052  →  **Cat.2 Refinamiento**
- **Estructura ya publicada:** El modelo Boltzmann–Vlasov con fuente de Schwinger y factor (1−2f) ya está en Kluger et al. (1991–98) y Prakapenia–Vereshchagin (2023); la frecuencia de plasma ∝ e√(eE) es la escala estándar.
- **Ecuación equivalente:** Ecuación de Vlasov cuántica en el límite de Markov con bloqueo de Pauli.
- **Aporte del repositorio:** Identidad simbólica de energía con (1−2f) en j_pol (K1) y solución exacta del bloqueo (K2); DOP853 con eventos; |ΔE|/E₀≤2.9e-06 tras ≥10 cruces; 0≤f≤½; constantes numéricas de ESTE modelo (ancho de fuente σ=0.25√ε₀; sensibilidad ≤0.1% en τ₀). Sin bloqueo f llega a 6.53; el bloqueo alarga τ₀ de 4.5576 a 5.5684 (ε₀=100).
- **No presentar como:** «Hemos calculado el amortiguamiento de las oscilaciones de plasma de Schwinger» como resultado de QED; las constantes ω_p/(e√ε₀) y γ/ω_p como universales.
- **Cómo citar:** «En un modelo de Vlasov de Markov con bloqueo de Pauli (Kluger et al.; Prakapenia–Vereshchagin) obtenemos ω_p≈0.19e√(eE₀) y γ/ω_p≈5·10⁻³ para una fuente gaussiana de ancho 0.25√ε₀.»
- **Salvedades:** El bloqueo de Markov acota f<½ desde el vacío, mientras que la ocupación exacta de campo constante e^{−πm_⊥²/ε} supera ½ para ε>4.5: el modelo subestima la ocupación a campo fuerte. Sin QKE no markoviana, sin transporte espacial, sin comparar con las curvas publicadas.

### P1 — Diagrama de fases (N, r₀)  →  **Cat.3 Síntesis**
- **Estructura ya publicada:** Cada frontera es estándar: r_s=2GM (Schwarzschild), gas de Fermi degenerado, ε=1 (Schwinger), densidad de Planck, N_crit=M_Pl²/m² (cota tipo Gibbons), Nα≳1 (supercritico, Zel'dovich–Greiner).
- **Ecuación equivalente:** —
- **Aporte del repositorio:** Colocarlas en un solo plano con valores explícitos, tiempos y verificación simbólica de cada umbral (P1–P9).
- **No presentar como:** «Diagrama de fases fundamental de la materia cargada»; la fase (ii) usa electrones degenerados que NO están ligados (Coulomb>Fermi para N≳6·10³).
- **Cómo citar:** «Mapa de escalas cuantitativo (herramienta organizativa), no una predicción de formación de agujeros negros.»
- **Salvedades:** Criterio de atrapamiento = condición de datos iniciales; la carrera dinámica no se resolvió espacialmente.

### T1 — Umbrales N*=2.74e+05, N_crit=5.71e+44, N_valid=9/α  →  **Cat.1 Redescubrimiento**
- **Estructura ya publicada:** N_crit: reformulación de la cota de masa de Gibbons para BH extremales; N_valid: supercriticidad de Z≈137–173. N* (τ_dis=r₀) depende de κ.
- **Ecuación equivalente:** —
- **Aporte del repositorio:** N* como cantidad de modelo (Cat.3 menor).
- **No presentar como:** «N* es un umbral físico universal».
- **Cómo citar:** «Orden de magnitud dependiente de modelo.»
- **Salvedades:** Ver κ.

### M1 — Metodología: pipeline «LLM propone, CAS ejecuta», gates simbólicos, controles negativos, informe autogenerado  →  **Cat.3 Síntesis**
- **Estructura ya publicada:** No es un resultado de física; valor de reproducibilidad.
- **Ecuación equivalente:** —
- **Aporte del repositorio:** 67 gates (15 controles negativos), trazabilidad número→clave JSON.
- **No presentar como:** Presentarlo como evidencia de corrección física.
- **Cómo citar:** «Software abierto de verificación».
- **Salvedades:** Un gate PASS verifica coherencia interna, no validez física (p. ej. el error del «−8.33e-4» pasó los gates).

## 4. Afirmaciones que NO pueden presentarse como «descubrimiento nuevo» (riesgo de rechazo)
- ~~«Descubrimos que la energía electrostática no cuenta para atrapar».~~ → «Recordamos que, en términos de la masa de Misner–Sharp M(R)=M−Q²/(2R), el criterio de atrapamiento es 2GM(R)≥R, de modo que la energía electrostática exterior no contribuye [Misner–Sharp 1964; Hod; Boulware 1973].»
- ~~«Ley nueva τ∝E^{-1/2}».~~ → «Consistente con la escala de plasma conocida, ω_p∼e√(eE) [Kluger et al.].»
- ~~«κ=13.5394 es una constante universal de la descarga de Schwinger».~~ → «En un modelo de Markov con LCFA y primer cero del campo como definición, obtenemos κ≈13.54.»
- ~~«El crecimiento logarítmico y el polo en log₁₀ε≈562 son resultados nuevos».~~ → «Reproducimos el límite logarítmico conocido del Lagrangiano de Heisenberg–Euler [Ritus; Lebedev–Ritus; Dunne].»
- ~~«Violación cuántica acotada de la DEC» como hallazgo; el valor −8.33e-4 como cota.~~ → «La parte real del Lagrangiano a un bucle viola la DEC en configuración eléctrica pura para ε>1.88; el margen relativo tiende a −α/3π, reflejo de la anomalía de traza.»
- ~~«La QED no censura el horizonte» como predicción física en r_Q.~~ → «Mostramos que la extrapolación cuártica del EEH produce un perfil cualitativamente distinto del de la integral de tiempo propio exacta; comparamos ambos en el sistema acoplado.»
- ~~«Nuevo modelo de retroacción».~~ → «Reimplementamos el modelo de Markov de Kluger–Mottola–Eisenberg con verificación simbólica del balance de energía.»
- ~~«Hemos calculado el amortiguamiento de las oscilaciones de plasma de Schwinger» como resultado de QED; las constantes ω_p/(e√ε₀) y γ/ω_p como universales.~~ → «En un modelo de Vlasov de Markov con bloqueo de Pauli (Kluger et al.; Prakapenia–Vereshchagin) obtenemos ω_p≈0.19e√(eE₀) y γ/ω_p≈5·10⁻³ para una fuente gaussiana de ancho 0.25√ε₀.»
- ~~«N* es un umbral físico universal».~~ → «Orden de magnitud dependiente de modelo.»

## 5. Debilidades que un árbitro señalaría (autocrítica)
1. **Ritus–Narozhny.** El informe afirma que el polo de Landau es inaccesible y que WEC/NEC/SEC valen hasta ε=10⁴⁰ con un bucle. Si la conjetura aplica (αχ^{2/3}≳1, χ≈1604), los lazos superiores son O(1) mucho antes. Debe presentarse como «a un bucle» y citar Fedotov (arXiv:1608.02261) y Mironov–Meuren–Fedotov (2020; identificador arXiv a confirmar). Aplicabilidad a E puro: sin verificar.
2. **Expansión derivativa / Drummond–Hathrell.** 04 resuelve el perfil hasta r_Q≈10⁻¹⁸ λ_C (ℓ≈1.8e-18 λ_C), donde la LCFA y el acoplo mínimo no están controlados (términos RF² de orden α(λ_C/L)²). Son soluciones formales del modelo, no QED.
3. **Parte imaginaria y estabilidad.** Para ε>1 la tasa de decaimiento del vacío es O(1) en unidades Compton; la parte real aislada no define una acción efectiva estable. El perfil estático para ε≫1 es un ejercicio de modelo.
4. **Causalidad (Russo–Townsend).** Causalidad⇒DEC y SEC en clase Plebanski; nuestra DEC se viola para ε>1.88. Hay que calcular las condiciones de causalidad con L₁(S,P) y ver si el efecto es una propiedad real o del recorte P=0.
5. **Hipótesis de Penrose cargada.** Ahora DERIVADA internamente (P10), pero sólo para datos simétricos en el tiempo, ρ≥0 y sin materia cargada exterior (como KWY). Malec–Ó Murchadha avisan de rupturas en colapso temprano; la dinámica con carga que sale NO está cubierta.
6. **Cinética con bloqueo de Pauli.** La ecuación de Markov con (1−2f) acota f<½ desde el vacío, mientras que la QKE exacta de campo constante da f=e^{−πm⊥²/ε}>½ a ε>4.5; los resultados de 05b (ω_p, γ) son de modelo.
7. **Modelo de materia de la fase (ii).** El gas de electrones degenerados no está ligado (Coulomb 4.2e+42× gravedad); r_s es una condición de datos iniciales.
8. **Constante κ.** Dependiente del modelo (n=1, sin (1−2f), LCFA, primer cero). Sin comparación con las curvas de KME/Bloch/Ruffini en texto completo.
9. **Verificación bibliográfica.** Ningún texto completo leído (arxiv.org y ar5iv no accesibles). Todos los niveles E son extractos; los M no se comprobaron.
10. **Error propio detectado.** DEC: el «−8.33e-4 acotado» era artefacto del rango (ε=10⁴⁰); umbral 1.92 era artefacto de rejilla (exacto 1.88218). Corregido en 04, informe y JSON.

## 6. El «hook» honesto
**Enunciado.** Un mapa de validez cuantitativo y verificable del sistema Einstein–Euler–Heisenberg: dónde la truncación cuártica es cualitativamente engañosa, dónde hace falta el Lagrangiano de un bucle exacto (tiempo propio), y dónde la EFT deja de existir (Schwinger ε~1, expansión derivativa/Compton, Ritus–Narozhny, Planck), junto con el contraste entre las escalas de descarga y de colapso.

**Por qué es honesto.** Cada ingrediente es conocido; lo inédito (hasta donde llegó una búsqueda superficial) es la comparación cuantitativa directa truncado-vs-exacto dentro del sistema acoplado, con las fronteras de validez en un mismo plano y con verificación formal reproducible.

**Contribuciones candidatas (todas condicionadas al §7):**
- La extrapolación cuártica del EEH a ε≫1 da E∼r^{-2/3}, masa interior finita y horizonte único tipo Schwarzschild; el tiempo propio exacto da RN con carga corrida (δr₊/r₊ medido vs predicho con 3% de acuerdo).
- Umbral exacto ε=1.88218 de inicio de violación de la DEC en configuración eléctrica pura y su identificación con la anomalía de traza (identidad S10).
- Mapa de escalas (N, r₀) con los tiempos τ_dis y τ_ff y sus cocientes, declarando explícitamente los supuestos.
- Software de verificación formal con controles negativos.

**Dónde.** Nota corta en gr-qc/hep-th (arXiv) como «validity map»; o artículo metodológico/de software (física computacional). NO como artículo de «descubrimiento» de censura cuántica.

**Lo que NO se afirma:** Que la QED censure o impida el colapso gravitatorio (ya implicado por Misner–Sharp/Boulware/Gibbons).; Leyes nuevas de plasma o de campo fuerte.; Predicciones físicas en r≪λ_C..

## 7. Verificaciones pendientes antes de cualquier envío
- [ ] Leer en texto completo: Kluger et al. (1991, 1992, 1998), Ruffini–Vereshchagin–Xue §8, Gelis–Tanji, Ruffini–Wu–Xue, Yajima–Tamaki, Boulware, Malec–Ó Murchadha, KWY, Russo–Townsend, Fedotov.
- [ ] Comparar κ reproduciendo una curva publicada (misma definición de tiempo/longitud de oscilación).
- [ ] HECHO en 05b: factor (1−2f) y oscilaciones multi-cruce. PENDIENTE: ecuación cinética cuántica no markoviana (f hasta 1) y comparación de ω_p, γ con Prakapenia–Vereshchagin (2023).
- [ ] Calcular las condiciones de causalidad de Russo–Townsend con L₁(S,P).
- [ ] Buscar sistemáticamente RN/EEH con carga corrida (Uehling) y logarítmico; no se hizo.
- [ ] Verificar contra la edición original los datos bibliográficos marcados VERIFY en references.bib (revista/volumen/páginas recordados de memoria) y añadir la afiliación en paper.tex (el autor ya figura: Álvaro López Almeida).
- [ ] Incluir los términos de Drummond–Hathrell y restringir 04 a r≳λ_C o rotular como modelo.
