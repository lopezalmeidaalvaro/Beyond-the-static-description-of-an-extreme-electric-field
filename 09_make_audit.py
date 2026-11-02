"""Genera NOVELTY_AUDIT.md y novelty_audit.json desde UNA estructura de datos + discovery_results.json.
Niveles de verificacion bibliografica (honestos):
  E = visto en extractos de busqueda esta sesion (resumenes/snippets; NO texto completo)
  M = de memoria del autor; NO verificado esta sesion
  N = no encontrado en las busquedas realizadas (ausencia en una busqueda superficial != inexistencia)
  T = texto completo leido -> NINGUNO (arxiv.org / ar5iv no accesibles desde esta sesion)"""
import json, datetime
D = json.load(open('discovery_results.json')); A, B, C, K = D['branch_A'], D['branch_B'], D['branch_C'], D['audit_checks']
ec = A['energy_conditions_qed']
V = dict(
  kappa=B['scaling']['kappa_inf'], r0s=C['race_law']['r0_star_over_lPlanck'], dec_on=ec['DEC_onset_eps_exact'],
  dec_asym=ec['DEC_margin_asymptote_alpha_over_3pi'], dec_1e10=ec['DEC_margin_at']['1e+10'], dec_1e40=ec['DEC_margin_at_range_edge_1e40'],
  landau_GeV=K['landau_scale_GeV_log10'], rn=K['ritus_narozhny_chi_threshold'], rn_orders=K['orders_above_RN_threshold_in_05'],
  ell=K['ell_over_lambdaC_at_rQ_N1e+06'], rr_ours10=K['ruffini_crosscheck']['ours']['10']['tau0_first_zero'][1],
  rr_ours015=K['ruffini_crosscheck']['ours']['0.15']['tau0_first_zero'][1], rr_ratio=K['ruffini_ratio_ours'],
  nstar=C['critical']['N_star'], ncrit=C['critical']['N_crit'], z2=C['critical']['z2_coulomb_over_gravity'],
  shift_m=A['QED_horizon_shift']['measured'], shift_p=A['QED_horizon_shift']['predicted_running_charge'],
  ngates=D['meta']['n_gates'], nneg=D['meta']['n_negative_controls'],
  om_ratio=D['branch_B2']['scaling']['omega_over_esqrtx0'][-1], gam_ratio=D['branch_B2']['scaling']['gamma_over_omega'][-1],
  emax=max(v['Emax_err'] for v in D['branch_B2']['sweep'].values()), tau_blk=D['branch_B2']['pauli_effect_x100']['tau0_block'],
  tau_nob=D['branch_B2']['pauli_effect_x100']['tau0_noblock'], fmax_nob=D['branch_B2']['pauli_effect_x100']['f_max_noblock'],
  dU=D['branch_C']['U_QED']['table']['1e+00']['delta_U'],
  b3_fmax=max(v['fmax'] for k, v in D['branch_B3']['exact_sweep'].items() if float(k) >= 10), b3_om=D['branch_B3']['scaling']['omega_over_esqrtx0'][-1],
  b3_gam=D['branch_B3']['scaling']['gamma_over_omega'][-1], b3_tau100=D['branch_B3']['exact_sweep']['100']['tau0_e_sqrtx0'], mk_tau100=D['branch_B3']['markov_sweep']['100']['tau0_e_sqrtx0'],
  b3_wpl=D['branch_B3']['exact_sweep']['100']['w_pl_avg_last_period'],
  c2_loss=D['branch_C2']['weak_gravity_F']['inst']['interior_only'], c2_gain=D['branch_C2']['weak_gravity_F']['inst']['exterior_only'], c2_net=D['branch_C2']['weak_gravity_F']['inst']['both'],
  c2_fesc=D['branch_C2']['weak_gravity_F']['inst']['f_esc_interior'], c2_shift=D['branch_C2']['threshold_range_pct'], c2_NA=D['branch_C2']['N_A_staticcurve_enters_initial_horizon'])

# ----------------------------------------------------------------------------------------------------------------
PAPERS = [
 dict(id='HE1936', ref='Heisenberg & Euler, Z. Phys. 98, 714 (1936); Weisskopf, Kgl. Danske Vid. Selsk. Mat.-fys. Medd. 14 (6) (1936)', lvl='M',
      know='Lagrangiano efectivo de QED a un bucle; término cuártico (2α²/45m⁴)[(E²−B²)²+7(E·B)²].', gates='04:S5,S6 (coef. 4cE⁴), 01/03 (baseline cuártico)', cat=1),
 dict(id='S1951', ref='Schwinger, Phys. Rev. 82, 664 (1951)', lvl='M',
      know='Forma de tiempo propio exacta; Im L = tasa de decaimiento del vacío; serie en e^{−nπ/ε}.', gates='04:N1 (Im L₁ = serie, 8e-15), S6; 05:D3 (tasa de pares)', cat=1),
 dict(id='R1975', ref='Ritus, Zh. Eksp. Teor. Fiz. 69, 1517 (1975); Ritus hep-th/9812124 (1998); Lebedev & Ritus (representación dispersiva)', lvl='E',
      know='Aparece como referencia estándar del Lagrangiano de campo intenso y QED a distancias cortas; un extracto de 2024 indica que el límite de campo fuerte es logarítmico con el coeficiente de la función β (Lebedev–Ritus: representación dispersiva). No se vio el texto.',
      gates='04:S7 (coef. −α/3π), S8 (polo), N3 (pendiente −1/24π²), N8', cat=1),
 dict(id='G1975', ref='Gibbons, Commun. Math. Phys. 44, 245 (1975)', lvl='E',
      know='Pérdida espontánea de carga de agujeros negros por emisión de partículas; para agujeros grandes (>~10¹⁷ g) rige una fórmula tipo Schwinger. Un trabajo de 2025 (arXiv:2512.01898) resume que Gibbons calculó la probabilidad de descarga con coeficientes de Bogoliubov para un campo escalar cargado, que los agujeros pequeños (campo local mayor) se descargan más, y que esto implica una cota inferior de masa para agujeros extremales.',
      gates='06:P5 (N_crit=M_Pl²/m²), 04 (ε(r_Q)=M_Pl²/N), §7 «facilitador»', cat=1),
 dict(id='DH1980', ref='Drummond & Hathrell, Phys. Rev. D 22, 343 (1980)', lvl='E',
      know='Según trabajos posteriores: acción efectiva a un bucle con términos RF², R_μνF², R_μνλρF² y un término derivativo, con coeficientes O(α)/m²; desplazamiento de velocidad de fotones O(α(λ_C/L)²), nulo en horizontes.',
      gates='NINGUNO: nuestro tratamiento curvo es de acoplo mínimo. LAGUNA: a L≈r_Q≪λ_C el cociente λ_C²/L² es enorme y la expansión no está controlada', cat=None),
 dict(id='YT2001', ref='Yajima & Tamaki, Phys. Rev. D 63, 064007 (2001), gr-qc/0005016', lvl='E',
      know='Agujeros negros estáticos esféricos en Einstein–Euler–Heisenberg (EH como teoría efectiva), caso magnético resuelto analíticamente; parámetros a,b; dependencia de causalidad y termodinámica; singularidad central; a>0 natural si EH es límite de baja energía de Born–Infeld.',
      gates='01/03 (EEH cuártico estático), 04:S5 (relación BI→EH)', cat=1),
 dict(id='RWX2013', ref='Ruffini, Wu & Xue, arXiv:1307.4951 (2013) [revista: Phys. Rev. D 88, 085004 – M]', lvl='E',
      know='Einstein + Lagrangiano EH a un bucle; casos eléctrico, magnético y diónico; horizonte, área, entropía, energía total y extraíble con corrección QED; la condición de agujero extremo se desplaza. Resultados a orden dominante en la corrección QED; el desarrollo O(α³) se escribe. Una búsqueda de seguimiento no encontró cálculo no perturbativo (tiempo propio exacto) de horizonte RN.',
      gates='01 (derivación), 03/04 (cuártico); 04 extiende a tiempo propio exacto (Cat.2)', cat=1),
 dict(id='CM1989', ref='Cooper & Mottola, Phys. Rev. D 40, 456 (1989) [volumen/página: M]', lvl='M',
      know='Retroacción cuántica como problema de valores iniciales (QED escalar). No localizado en las búsquedas.', gates='05 (contexto)', cat=1),
 dict(id='KESCM1991', ref='Kluger, Eisenberg, Svetitsky, Cooper, Mottola, PRL 67, 2427 (1991); PRD 45, 4659 (1992) [M]; Int. J. Mod. Phys. E (1993) = hep-ph/0311293', lvl='E',
      know='Producción de pares en campo eléctrico fuerte como problema de valores iniciales en campo medio semiclásico homogéneo; la retroacción genera una corriente que se opone al campo y aparecen oscilaciones de plasma. Se compara con un modelo Boltzmann–Vlasov con fuente de Schwinger y factor de bloqueo de Pauli (1−2f).',
      gates='05 (modelo completo: fuente, corriente de conducción + polarización, D2 balance de energía, oscilaciones)', cat=1),
 dict(id='KME1998', ref='Kluger, Mottola & Eisenberg, hep-ph/9803372 (1998) [PRD 58, 125015 – M]', lvl='E',
      know='El número de partículas adiabático obedece una ecuación de Vlasov cuántica no local en el tiempo; para campos débiles y lentos actúa como función de distribución y su derivada como fuente de Boltzmann–Vlasov (límite de Markov).',
      gates='05 (fuente local = límite de Markov; nosotros: solo término n=1, sin (1−2f), válido hasta el primer cero)', cat=1),
 dict(id='SCH1998', ref='Schmidt et al., hep-ph/9809227, Int. J. Mod. Phys. E 7 (1998)', lvl='E',
      know='Ecuación cinética cuántica no markoviana; a baja densidad el término fuente se reduce a la parte principal de la fórmula de Schwinger.', gates='05 (justificación de la fuente local)', cat=1),
 dict(id='RVX2010', ref='Ruffini, Vereshchagin & Xue, Phys. Rep. 487, 1 (2010), arXiv:0910.0974', lvl='E',
      know='Revisión; §8 trata oscilaciones de plasma, retroacción y apantallamiento por pares. Un extracto (de segunda mano) da longitudes de oscilación: E₀=10E_c → ~10 ħ/m_ec; E₀=0.15E_c → ~10⁷ ħ/m_ec.',
      gates=f"05: contraste cuantitativo (08): nuestro primer cero {V['rr_ours10']:.3g} y {V['rr_ours015']:.3g} (cociente {V['rr_ratio']:.2g} vs 10⁶ citado)", cat=1),
 dict(id='GT2016', ref='Gelis & Tanji, Prog. Part. Nucl. Phys. 87, 1 (2016), arXiv:1510.05451', lvl='E',
      know='Solo metadatos y resumen (revisión del mecanismo de Schwinger, incl. Schwinger asistido dinámicamente). Contenido NO accedido.', gates='05 (no comparado; pendiente)', cat=None),
 dict(id='PV2023', ref='Prakapenia & Vereshchagin, arXiv:2306.12956 (2023)', lvl='E',
      know='El bloqueo de Pauli suprime la creación de pares en campos sobrecríticos y modifica la retroacción; cuando los pares dominan energéticamente la frecuencia de oscilación aumenta.', gates='05: límite declarado (validez hasta el primer cero)', cat=1),
 dict(id='B1973', ref='Boulware, Phys. Rev. D 8, 2363 (1973)', lvl='E',
      know='Fuente secundaria: capa fina cargada con uniones de Israel y extensión máxima de RN; ecuación de energía que relaciona masa total, masa propia, carga y radio; singularidad desnuda solo con masa propia negativa.', gates='06:P1–P3 (L1), §7', cat=1),
 dict(id='BEK1971', ref='Bekenstein, Phys. Rev. D 4, 2185 (1971)', lvl='E',
      know='Equilibrio hidrostático y colapso de bolas de fluido cargado; se cita como origen del «doble papel» de la carga y del teorema de Birkhoff cargado. Contenido detallado NO accedido.', gates='06 (L1, z²), §7', cat=1),
 dict(id='MOM1994', ref='Malec & Ó Murchadha, Phys. Rev. D 49 (1994) [M], gr-qc/9401014', lvl='E',
      know='Desigualdad de Penrose en simetría esférica aun con el horizonte dentro de la materia; condición necesaria para formación de superficies atrapadas; la versión cargada tipo Gibbons se demuestra si toda la carga está dentro de la superficie atrapada exterior y, según el resumen, puede romperse en etapas tempranas del colapso.', gates='06:P1–P3, criterio (ii); riesgo: hipótesis de carga contenida', cat=1),
 dict(id='KWY', ref='Khuri, Weinstein & Yamada: arXiv:1308.3771 (2013), arXiv:1409.3271 (2014/15), Contemp. Math. (2015)', lvl='E',
      know='Desigualdad de Penrose riemanniana con carga (cota superior r ≤ m+√(m²−q²)) bajo DEC cargada y SIN materia cargada fuera del horizonte; m≥|q|. El «2017» del encargo corresponde, de memoria (NO verificado), a la versión publicada de arXiv:1409.3271 en J. Differential Geom. 106 (2017); lo visto en extractos es 2013–2015 (preprints/Contemp. Math.).', gates='06:P1 (RN satura), criterio (ii)', cat=1),
 dict(id='MS1964', ref='Misner & Sharp, Phys. Rev. 136, B571 (1964) [M]; Hayward (1994) trapping horizons [E parcial]; Hod (masa M∞−Q²/2R) [E vía arXiv:2606.31116]', lvl='E',
      know='Un trabajo reciente indica que, en RN, la masa dentro de R, M(r≤R)=M∞−Q²/2R, coincide exactamente con la masa cuasilocal de Misner–Sharp M=(R/2)(1−g^{rr}); el horizonte de atrapamiento es el borde 2M(R)=R.', gates='06:P2, P3 = ESTA identidad. L1 es una reescritura.', cat=1),
 dict(id='RT2024', ref='Russo & Townsend, JHEP 06 (2024) 191, arXiv:2404.09994', lvl='E',
      know='En NLED de clase Plebanski, causalidad implica DEC y SEC; toda violación implica acausalidad. Tomizawa & Suzuki (arXiv:2309.10535): violación de DEC ↔ fotones superlumínicos en la geometría efectiva.', gates='04:N8, N9, S10 (nuestra DEC solo en configuración eléctrica pura, P=0): TENSIÓN a comprobar', cat=None),
 dict(id='RNC', ref='Fedotov et al., arXiv:1608.02261; Mironov, Meuren & Fedotov (2020) [identificador arXiv:2003.06909 NO confirmado]; Heinzl, Ilderton & King (unitariedad)', lvl='E',
      know=f"Conjetura de Ritus–Narozhny: cada lazo aporta αχ^{{2/3}}, no α; la teoría perturbativa fallaría para αχ^{{2/3}}≳1 (χ≳α^(-3/2)≈{V['rn']:.0f}) en campo cruzado constante. La resumación indica ruptura ya en g≳1; Heinzl–Ilderton–King discuten que la unitariedad elimina diagramas supuestamente responsables. Abierta; se enuncia para campo cruzado constante (aplicabilidad a E puro: NO verificada).",
      gates='04:N7,N8; 05 barrido hasta ε=10¹⁰; DISCOVERY_REPORT §4/§6. LAGUNA importante', cat=None),
 dict(id='ZGP', ref='Gershtein & Zel\'dovich (1969); Pieper & Greiner (1969); Greiner, Müller & Rafelski (1985) [M]; Popov et al. arXiv:2008.05005, 1910.01373 [E]', lvl='E',
      know='Núcleo extendido: el nivel 1s alcanza el continuo negativo en Z≈173; carga puntual: el nivel existe solo hasta Z≈137; decaimiento del vacío con emisión espontánea de positrones; apantallamiento y bloqueo de Pauli reducen la producción.', gates='06: N_valid=9/α, «r_neutral», supercriticidad Nα≳1', cat=1),
 dict(id='TA1977', ref='Adler, Collins & Duncan, Phys. Rev. D 15, 1712 (1977); Collins, Duncan & Joglekar, Phys. Rev. D 16, 438 (1977) [M]', lvl='M',
      know=f"Anomalía de traza T^μ_μ=(β/2e)F². En nuestro resultado: margen DEC asintótico −α/(3π)={V['dec_asym']:.3e}.", gates='04:S10, N9', cat=1),
 dict(id='XBH', ref='Worldline-instanton / Schwinger en RN extremal: arXiv:2409.02197; revisita de Gibbons: arXiv:2512.01898; arXiv:1712.01302', lvl='E',
      know='La tasa de Schwinger en RN extremal se apaga al acercarse la razón q/m a la cota de extremalidad; el efecto es máximo en el horizonte y decae en una escala Compton; un RN se descarga apreciablemente cuando m² ≲ (carga·campo)/(radio transversal)².', gates='06 (N_crit), §7', cat=1),
]

RESULTS = [
 dict(id='L1', name='Cancelación de Coulomb en el atrapamiento de un cascarón (P1–P3)', cat=1,
      evidence='Misner–Sharp en RN: M(R)=M∞−Q²/2R (Hod; Misner & Sharp 1964); ecuación de energía de capa fina (Israel; Boulware 1973); desigualdad de Penrose cargada saturada por RN (Gibbons 1972 [M]; Malec–Ó Murchadha 1994; Khuri–Weinstein–Yamada).',
      equivalent='G·(M−Q²/8πR)=R/2 ⇔ g^{rr}(R)=0 (definición del horizonte de RN).',
      ours='CIERRE (P10, P11): derivación interna con SymPy de θ₊=0 ⇔ 2GM_MS=R, M_ADM=M_MS+U_C y la igualdad de Penrose cargada para datos simétricos en el tiempo, sin teoremas externos; U_C exacta a un bucle (δU/U_C=%.3e en R=r_Q, N=1e6). Sigue siendo Cat.1: el contenido es conocido; lo que cambia es que ya no depende de una cita.' % V['dU'],
      forbidden='«Descubrimos que la energía electrostática no cuenta para atrapar».',
      wording='«Recordamos que, en términos de la masa de Misner–Sharp M(R)=M−Q²/(2R), el criterio de atrapamiento es 2GM(R)≥R, de modo que la energía electrostática exterior no contribuye [Misner–Sharp 1964; Hod; Boulware 1973].»',
      caveat='Probado internamente sólo para datos esféricos simétricos en el tiempo, ρ≥0 y sin materia cargada exterior (las mismas hipótesis de KWY); el caso dinámico general (Misner–Sharp completo) NO. Malec–Ó Murchadha indican que la versión cargada puede romperse en etapas tempranas del colapso si hay carga exterior.'),
 dict(id='L2a', name='Escala τ_dis ∝ (eE₀)^{-1/2} (exponente −0.5000)', cat=1,
      evidence='Análisis dimensional con m→0 (única escala √(eE)); frecuencia de plasma ω_p∝e√(eE) en la retroacción de Schwinger (Kluger et al. 1991–1998; Ruffini–Vereshchagin–Xue 2010 §8) [fórmula exacta: verificar en texto completo].',
      equivalent='t·√(eE₀)=κ(e) por invariancia de escala.', ours='Medición numérica del exponente y verificación de convergencia.',
      forbidden='«Ley nueva τ∝E^{-1/2}».', wording='«Consistente con la escala de plasma conocida, ω_p∼e√(eE) [Kluger et al.].»', caveat='Exponente = consecuencia dimensional, no resultado.'),
 dict(id='L2b', name=f"Constante κ = {V['kappa']:.4f}", cat=2,
      evidence=f"Constante numérica de un modelo concreto (fuente local n=1, LCFA, primer cero del campo, 8 nodos Gauss–Laguerre en p⊥). Contraste con segunda mano: primeros ceros {V['rr_ours10']:.3g} (E₀=10E_c) y {V['rr_ours015']:.3g} (0.15E_c) frente a ≈10 y ≈10⁷ en un extracto de Ruffini et al.; coinciden en orden de magnitud en ambos extremos (cociente {V['rr_ratio']:.2g} vs 10⁶).",
      equivalent='—', ours='Valor numérico con convergencia (orden 1.98) y método independiente.',
      forbidden='«κ=13.5394 es una constante universal de la descarga de Schwinger».',
      wording='«En un modelo de Markov con LCFA y primer cero del campo como definición, obtenemos κ≈13.54.»',
      caveat='Dependencia de modelo O(1); sin (1−2f) ni QKE; definición «primer cero» no estándar; comparación incompleta con la literatura.'),
 dict(id='L3', name=f"Ley τ_dis/τ_ff=√(r₀*/r₀), r₀*={V['r0s']:.2e} l_P, independiente de N", cat=3,
      evidence='Es el producto de dos leyes de potencia con el mismo √N (τ_dis∝r₀/√N, τ_ff∝√(r₀³/N)); r₀* está fijada por e²/(Gm²)=αM_Pl² (el cociente EM/gravedad del electrón, conocido).',
      equivalent=f"τ_dis/τ_ff ∼ (Gm²/e²)^{{1/2}}·(κ-factores); z²={V['z2']:.2e}.", ours='Prefactor cerrado con κ y verificación simbólica (P7–P9).',
      forbidden='«Ley universal sin precedentes»; r₀*≪l_P está fuera de la validez de ambos modelos (solo la desigualdad τ_dis<τ_ff para r₀>l_P tiene sentido).',
      wording='«Para N electrones, el cociente entre el tiempo de descarga y el de caída libre escala como (r₀*/r₀)^{1/2} con r₀* determinado por αM_Pl² y κ.»',
      caveat='A r₀≲λ_C la caída libre de polvo de N masas en reposo no describe materia degenerada confinada.'),
 dict(id='L4a', name='Crecimiento logarítmico a un bucle y polo de Landau', cat=1,
      evidence=f"Weisskopf 1936; Schwinger 1951; Ritus 1975; Lebedev–Ritus; coeficiente de la función β. Cross-check cuantitativo: nuestro ε_L → μ_L≈10^{V['landau_GeV']:.1f} GeV frente a ≈10²⁷⁷ GeV citado para el polo de Landau de la QED (extracto de arXiv:1812.08673).",
      equivalent='L=(E²/2)[1−(α/3π)(ln ε+const)].', ours='Constante asintótica K y evaluación numérica exacta del rayo (Re y Im).',
      forbidden='«El crecimiento logarítmico y el polo en log₁₀ε≈562 son resultados nuevos».', wording='«Reproducimos el límite logarítmico conocido del Lagrangiano de Heisenberg–Euler [Ritus; Lebedev–Ritus; Dunne].»',
      caveat=f"Solo un bucle. Conjetura de Ritus–Narozhny: αχ^{{2/3}}≳1 ya en χ≈{V['rn']:.0f}; nuestro barrido llega a ε=10¹⁰ ({V['rn_orders']:.1f} órdenes por encima)."),
 dict(id='L4b', name=f"Violación de DEC desde ε={V['dec_on']:.5f}; margen {V['dec_asym']:.3e} (asintótico)", cat=1,
      evidence='Anomalía de traza a un bucle (Adler–Collins–Duncan 1977; Collins–Duncan–Joglekar 1977 [M]). Identidad S10: E·D−2L=x³d(L₁/x²)/dx.',
      equivalent=f"margen→−α/(3π)={V['dec_asym']:.3e}; a ε=10¹⁰: {V['dec_1e10']:.3e}.",
      ours=f"Umbral exacto ε={V['dec_on']:.5f} (máximo de L₁/x²) con integral exacta (Cat.2, menor). CORRECCIÓN: el «−8.33e-4 acotado» del informe previo era el borde del rango (ε=10⁴⁰: {V['dec_1e40']:.3e}); el margen NO está acotado.",
      forbidden='«Violación cuántica acotada de la DEC» como hallazgo; el valor −8.33e-4 como cota.',
      wording='«La parte real del Lagrangiano a un bucle viola la DEC en configuración eléctrica pura para ε>1.88; el margen relativo tiende a −α/3π, reflejo de la anomalía de traza.»',
      caveat='Tensión con Russo–Townsend (causalidad⇒DEC): falta calcular las condiciones de causalidad con L₁(S,P) completo; para ε>1 Im L≠0, la parte real aislada no es una acción efectiva estable.'),
 dict(id='A1', name='Perfil interior estático: EH cuártico vs un bucle exacto vs Born–Infeld (04)', cat=2,
      evidence='Yajima–Tamaki 2001 y Ruffini–Wu–Xue 2013 trabajan con EH truncado (cuártico/leading order); una búsqueda de seguimiento no halló horizonte de RN con el tiempo propio exacto en campo fuerte (nivel N).',
      equivalent='—', ours=f"E∼r^{{-2/3}} y masa finita del cuártico como artefacto (S9); un bucle exacto = RN con carga corrida (δr₊/r₊ medido {V['shift_m']:.3e} vs predicho {V['shift_p']:.3e}); BI/cuártico → horizonte único tipo Schwarzschild.",
      forbidden='«La QED no censura el horizonte» como predicción física en r_Q.', wording='«Mostramos que la extrapolación cuártica del EEH produce un perfil cualitativamente distinto del de la integral de tiempo propio exacta; comparamos ambos en el sistema acoplado.»',
      caveat=f"En los perfiles de 04, la escala de variación del campo es ℓ≈r_Q/2={V['ell']:.1e} λ_C ≪ λ_C: viola el criterio derivativo del propio baseline (03) y la expansión de Drummond–Hathrell no está controlada. Son soluciones FORMALES del modelo de campo localmente constante."),
 dict(id='B1', name='Modelo cinético 05 (fuente, conducción + polarización, balance de energía)', cat=1,
      evidence='Kluger et al. 1991–98 (modelo y límite de Markov), Schmidt et al. 1998.', equivalent='Mismas ecuaciones esenciales.',
      ours='Verificación simbólica D1/D2, implementación Volterra con convergencia de orden 2 y controles negativos.', forbidden='«Nuevo modelo de retroacción».', wording='«Reimplementamos el modelo de Markov de Kluger–Mottola–Eisenberg con verificación simbólica del balance de energía.»', caveat='Sin (1−2f), sin QKE, sin transporte espacial; válido hasta el primer cero.'),
 dict(id='K1', name=f"Cinética de Markov con bloqueo de Pauli a través de ≥5 oscilaciones (05b): ω_p={V['om_ratio']:.4f} e√ε₀, γ/ω_p={V['gam_ratio']:.4f}", cat=2,
      evidence='El modelo Boltzmann–Vlasov con fuente de Schwinger y factor (1−2f) ya está en Kluger et al. (1991–98) y Prakapenia–Vereshchagin (2023); la frecuencia de plasma ∝ e√(eE) es la escala estándar.',
      equivalent='Ecuación de Vlasov cuántica en el límite de Markov con bloqueo de Pauli.',
      ours=f"Identidad simbólica de energía con (1−2f) en j_pol (K1) y solución exacta del bloqueo (K2); DOP853 con eventos; |ΔE|/E₀≤{V['emax']:.1e} tras ≥10 cruces; 0≤f≤½; constantes numéricas de ESTE modelo (ancho de fuente σ=0.25√ε₀; sensibilidad ≤0.1% en τ₀). Sin bloqueo f llega a {V['fmax_nob']:.2f}; el bloqueo alarga τ₀ de {V['tau_nob']:.4f} a {V['tau_blk']:.4f} (ε₀=100).",
      forbidden='«Hemos calculado el amortiguamiento de las oscilaciones de plasma de Schwinger» como resultado de QED; las constantes ω_p/(e√ε₀) y γ/ω_p como universales.',
      wording='«En un modelo de Vlasov de Markov con bloqueo de Pauli (Kluger et al.; Prakapenia–Vereshchagin) obtenemos ω_p≈0.19e√(eE₀) y γ/ω_p≈5·10⁻³ para una fuente gaussiana de ancho 0.25√ε₀.»',
      caveat='El bloqueo de Markov acota f<½ desde el vacío, mientras que la ocupación exacta de campo constante e^{−πm_⊥²/ε} supera ½ para ε>4.5: el modelo subestima la ocupación a campo fuerte. Sin QKE no markoviana, sin transporte espacial, sin comparar con las curvas publicadas.'),
 dict(id='E1', name=f"Cinética cuántica exacta (Bloch/Bogoliubov acoplada a Maxwell, 05c): f_max→{V['b3_fmax']:.4f}, ω_p={V['b3_om']:.4f} e√ε₀, γ/ω_p={V['b3_gam']:.4f}, ⟨w∥⟩_pl={V['b3_wpl']:.3f}", cat=1,
      evidence='Las ecuaciones de Bloch/Bogoliubov de la QKE para Dirac en E(t) homogéneo y su acoplo de Maxwell son estándar (Kluger et al. 1991–98; Schmidt et al. 1998; Prakapenia–Vereshchagin 2023 las integran con retroacción y bloqueo de Pauli exacto).',
      equivalent='QKE no markoviana en forma de Bloch con retroacción de Maxwell.',
      ours=f"Derivación simbólica Dirac→Bloch (Q0), pureza exacta (Q1), energía exacta con corriente de polarización (Q2), traza (Q4), Landau–Zener numérico vs e^(−πm⊥²/eE) (Q3), núcleo C verificado; barrido ε₀∈[4,10³] con ≥5 oscilaciones (ε₀≥10), convergencia en P, p⊥, Δt, rango y estado inicial; comparación directa con Markov (τ₀·e√ε₀: {V['b3_tau100']:.2f} vs {V['mk_tau100']:.2f} a ε₀=100). Detectado un factor 2 en el acoplo del encargo (C11).",
      forbidden='«Primera solución no markoviana de las oscilaciones de plasma de Schwinger»; «ley universal ω_p=0.19e√(eE₀)» (el valor exacto es otro y depende de la condición inicial y de la ausencia de colisiones).',
      wording='«Integramos las ecuaciones de Bloch de la QKE con retroacción (como Kluger et al. y Prakapenia–Vereshchagin) y cuantificamos la diferencia con el límite de Markov con bloqueo de Pauli en el mismo sistema.»',
      caveat='Modelo homogéneo, sin colisiones ni transporte (el plasma NO se isotropiza: w∥→≈1, no 1/3); estado inicial = vacío adiabático de orden cero (ambigüedad cuantificada, K8); no se ha reproducido ninguna curva publicada (texto completo no accesible).'),
 dict(id='DM1', name="Balance de Misner–Sharp dinámico (K_ij≠0) y transporte de flujos nulos radiales (06b: D_MS1–D_MS5, C_MS)", cat=1,
      evidence='Las ecuaciones ∂_rM=4πr²ρ, ∂_tM=4πr²T^r_t y el balance de flujo cuasilocal son estándar (Misner–Sharp 1964; Hayward 1994); el transporte de polvo nulo radial y su colapso (Vaidya, Wheeler) son clásicos.',
      equivalent='Ecuaciones de Misner–Sharp con flujo; polvo nulo en calibre polar-areal.',
      ours='Derivación completa en SymPy con φ(t,r), λ(t,r), campo de Maxwell y término de intercambio F_νλJ^λ; la conversión interna cancela y solo el flujo en r=R cambia M; control con ley de transferencia alterada; vínculo de momento preservado por el esquema; residuo exacto de la materia mantenida.',
      forbidden='«Primera demostración del balance de Misner–Sharp con descarga»; «resolvimos el colapso cargado dinámico completo».',
      wording='«Verificamos simbólicamente, en la métrica general esférica, que la conversión local campo→plasma conserva la masa de Misner–Sharp y que solo el flujo por la superficie la cambia.»',
      caveat='Es una verificación de identidades conocidas en un CAS, no un resultado nuevo de relatividad.'),
 dict(id='DM2', name=f"Descarga espacialmente resuelta de una bola cargada (06b): fuga {V['c2_loss']:+.3f} q², inflow {V['c2_gain']:+.3f} q², neto {V['c2_net']:+.3f} q²; umbral desplazado ≤{V['c2_shift']:.1f}%", cat=3,
      evidence='Cada ingrediente es conocido (descarga de Schwinger, colapso de radiación, extremalidad M≥Q de la bola cargada); la combinación cuantitativa con la curva exacta de 05c no se encontró en la búsqueda superficial.',
      equivalent='—',
      ours=f"Solver 1+1D validado (límite plano analítico, balance a 1e-3, convergencia, controles negativos); la fracción del campo interior que escapa por r₀ es {V['c2_fesc']:.2f}, pero el inflow del campo exterior es ~2.5 veces mayor; la dinámica en t≤r₀/c desplaza el umbral de horizonte <{V['c2_shift']:.1f}%; la curva r_s(N) del paper cae dentro del horizonte aparente inicial para N>{V['c2_NA']:.1e}.",
      forbidden='«La descarga de Schwinger acelera/permite el colapso» o «no lo retrasa» como predicción física; resultados para t>r₀/c.',
      wording='«En un modelo de flujos nulos radiales con cierre local de descarga, la dinámica de descarga hasta el tiempo de cruce de luz desplaza el umbral de horizonte de la bola cargada en menos de 1% hasta el tiempo de caída libre prescrito.»',
      caveat='Flujos radiales con núcleo reflectante (regularización de la cáustica), materia como masa de prueba y tiempo de caída libre PRESCRITO (el colapso autoconsistente de la materia cargada no se evoluciona), sin explosión de Coulomb ni colisiones, cierre de descarga local, sin Gauss durante la conversión; no se comparó con ningún cálculo publicado.'),
 dict(id='DM3', name="Colapso de un cascarón cargado con descarga: ley s_c = A q²/C^{3/2} y carrera contra el empuje de Coulomb (06b: P13, S1–S5)", cat=1,
      evidence='La unión de Israel y el colapso de cascarones cargados (Israel 1966; Boulware 1973; Bekenstein 1971) y la condición Q<M para el colapso son clásicos; la explosión coulombiana de materia cargada es bien conocida.',
      equivalent='Cascarón de Israel con carga: colapsa si Q<√(4πG)M.',
      ours='Ley analítica de la carrera s_c=2√I q²/(κC^{3/2}) con la curva exacta de descarga de 05c, verificada por simulación de la dinámica completa (error ~2%); evaluación en el q/C físico (≈1e21); la ganancia de masa de MS ΔE/E₀=(C/4)(s_c/s)²; controles negativos; integral primera conservada a 1e-12.',
      forbidden='«La QED/descarga de Schwinger impide (o permite) el colapso gravitatorio» como predicción física; N_race≈1e68 como resultado cuantitativo fiable (depende de un cierre local de la descarga y de la neutralización permanente).',
      wording='«En un modelo de cascarón fino con apantallamiento local de la descarga, el empuje de Coulomb desliga al cascarón salvo que el apantallamiento supere s_c=0.66 q²/C^{3/2}; con la relación carga-masa del electrón ese umbral excede con creces el disponible.»',
      caveat='Cierre local de la descarga y neutralización permanente tras el primer cero (ansatz); plasma y cascarón no acoplados; interior plano; una sola capa; sin explosión coulombiana de muchos cuerpos ni estructura radial; sin comparación con literatura (texto completo no accesible).'),
 dict(id='P1', name='Diagrama de fases (N, r₀)', cat=3,
      evidence='Cada frontera es estándar: r_s=2GM (Schwarzschild), gas de Fermi degenerado, ε=1 (Schwinger), densidad de Planck, N_crit=M_Pl²/m² (cota tipo Gibbons), Nα≳1 (supercritico, Zel\'dovich–Greiner).',
      equivalent='—', ours='Colocarlas en un solo plano con valores explícitos, tiempos y verificación simbólica de cada umbral (P1–P9).',
      forbidden='«Diagrama de fases fundamental de la materia cargada»; la fase (ii) usa electrones degenerados que NO están ligados (Coulomb>Fermi para N≳6·10³).', wording='«Mapa de escalas cuantitativo (herramienta organizativa), no una predicción de formación de agujeros negros.»',
      caveat='Criterio de atrapamiento = condición de datos iniciales; la carrera dinámica no se resolvió espacialmente.'),
 dict(id='T1', name=f"Umbrales N*={V['nstar']:.3g}, N_crit={V['ncrit']:.3g}, N_valid=9/α", cat=1,
      evidence='N_crit: reformulación de la cota de masa de Gibbons para BH extremales; N_valid: supercriticidad de Z≈137–173. N* (τ_dis=r₀) depende de κ.', equivalent='—', ours='N* como cantidad de modelo (Cat.3 menor).', forbidden='«N* es un umbral físico universal».', wording='«Orden de magnitud dependiente de modelo.»', caveat='Ver κ.'),
 dict(id='M1', name='Metodología: pipeline «LLM propone, CAS ejecuta», gates simbólicos, controles negativos, informe autogenerado', cat=3,
      evidence='No es un resultado de física; valor de reproducibilidad.', equivalent='—', ours=f"{V['ngates']} gates ({V['nneg']} controles negativos), trazabilidad número→clave JSON.", forbidden='Presentarlo como evidencia de corrección física.', wording='«Software abierto de verificación».', caveat='Un gate PASS verifica coherencia interna, no validez física (p. ej. el error del «−8.33e-4» pasó los gates).'),
]

WEAK = [
 ('Ritus–Narozhny', f"El informe afirma que el polo de Landau es inaccesible y que WEC/NEC/SEC valen hasta ε=10⁴⁰ con un bucle. Si la conjetura aplica (αχ^{{2/3}}≳1, χ≈{V['rn']:.0f}), los lazos superiores son O(1) mucho antes. Debe presentarse como «a un bucle» y citar Fedotov (arXiv:1608.02261) y Mironov–Meuren–Fedotov (2020; identificador arXiv a confirmar). Aplicabilidad a E puro: sin verificar."),
 ('Expansión derivativa / Drummond–Hathrell', f"04 resuelve el perfil hasta r_Q≈10⁻¹⁸ λ_C (ℓ≈{V['ell']:.1e} λ_C), donde la LCFA y el acoplo mínimo no están controlados (términos RF² de orden α(λ_C/L)²). Son soluciones formales del modelo, no QED."),
 ('Parte imaginaria y estabilidad', 'Para ε>1 la tasa de decaimiento del vacío es O(1) en unidades Compton; la parte real aislada no define una acción efectiva estable. El perfil estático para ε≫1 es un ejercicio de modelo.'),
 ('Causalidad (Russo–Townsend)', 'Causalidad⇒DEC y SEC en clase Plebanski; nuestra DEC se viola para ε>1.88. Hay que calcular las condiciones de causalidad con L₁(S,P) y ver si el efecto es una propiedad real o del recorte P=0.'),
 ('Hipótesis de Penrose cargada', 'Ahora DERIVADA internamente (P10), pero sólo para datos simétricos en el tiempo, ρ≥0 y sin materia cargada exterior (como KWY). Malec–Ó Murchadha avisan de rupturas en colapso temprano; la dinámica con carga que sale NO está cubierta.'),
 ('Cinética con bloqueo de Pauli', 'La ecuación de Markov con (1−2f) acota f<½ desde el vacío, mientras que la QKE exacta de campo constante da f=e^{−πm⊥²/ε}>½ a ε>4.5; los resultados de 05b (ω_p, γ) son de modelo. 05c resuelve esa limitación dentro del modelo homogéneo (QKE de Bloch exacta), cuantifica la diferencia, y NO incluye colisiones, transporte ni comparación con curvas publicadas.'),
 ('Brecha espacial/gravitatoria (06b)', 'El solver 1+1D llega a (π/2) r₀ solo gracias a un núcleo reflectante que regulariza la cáustica de los flujos radiales; la materia es masa de prueba y el tiempo de caída libre está prescrito: el colapso autoconsistente de un cascarón cargado (tensión radial, fuerza de Coulomb, salto de lapso) sigue sin resolverse; el desplazamiento <0.5% es del umbral de horizonte aparente, no de r_s del cascarón.'),
 ('Modelo de materia de la fase (ii)', f"El gas de electrones degenerados no está ligado (Coulomb {V['z2']:.1e}× gravedad); r_s es una condición de datos iniciales."),
 ('Constante κ', 'Dependiente del modelo (n=1, sin (1−2f), LCFA, primer cero). Sin comparación con las curvas de KME/Bloch/Ruffini en texto completo.'),
 ('Verificación bibliográfica', 'Ningún texto completo leído (arxiv.org y ar5iv no accesibles). Todos los niveles E son extractos; los M no se comprobaron.'),
 ('Error propio detectado', f"DEC: el «−8.33e-4 acotado» era artefacto del rango (ε=10⁴⁰); umbral 1.92 era artefacto de rejilla (exacto {V['dec_on']:.5f}). Corregido en 04, informe y JSON."),
]

HOOK = dict(
 statement='Un mapa de validez cuantitativo y verificable del sistema Einstein–Euler–Heisenberg: dónde la truncación cuártica es cualitativamente engañosa, dónde hace falta el Lagrangiano de un bucle exacto (tiempo propio), y dónde la EFT deja de existir (Schwinger ε~1, expansión derivativa/Compton, Ritus–Narozhny, Planck), junto con el contraste entre las escalas de descarga y de colapso.',
 why_honest='Cada ingrediente es conocido; lo inédito (hasta donde llegó una búsqueda superficial) es la comparación cuantitativa directa truncado-vs-exacto dentro del sistema acoplado, con las fronteras de validez en un mismo plano y con verificación formal reproducible.',
 candidate_claims=[
   'La extrapolación cuártica del EEH a ε≫1 da E∼r^{-2/3}, masa interior finita y horizonte único tipo Schwarzschild; el tiempo propio exacto da RN con carga corrida (δr₊/r₊ medido vs predicho con 3% de acuerdo).',
   f"Umbral exacto ε={V['dec_on']:.5f} de inicio de violación de la DEC en configuración eléctrica pura y su identificación con la anomalía de traza (identidad S10).",
   'Mapa de escalas (N, r₀) con los tiempos τ_dis y τ_ff y sus cocientes, declarando explícitamente los supuestos.',
   'Software de verificación formal con controles negativos.'],
 venues='Nota corta en gr-qc/hep-th (arXiv) como «validity map»; o artículo metodológico/de software (física computacional). NO como artículo de «descubrimiento» de censura cuántica.',
 not_claimed=['Que la QED censure o impida el colapso gravitatorio (ya implicado por Misner–Sharp/Boulware/Gibbons).',
              'Leyes nuevas de plasma o de campo fuerte.', 'Predicciones físicas en r≪λ_C.'],
 before_submission=['Leer en texto completo: Kluger et al. (1991, 1992, 1998), Ruffini–Vereshchagin–Xue §8, Gelis–Tanji, Ruffini–Wu–Xue, Yajima–Tamaki, Boulware, Malec–Ó Murchadha, KWY, Russo–Townsend, Fedotov.',
                    'Comparar κ reproduciendo una curva publicada (misma definición de tiempo/longitud de oscilación).',
                    'HECHO en 05b/05c: factor (1−2f), oscilaciones multi-cruce y QKE de Bloch exacta (f hasta 1). PENDIENTE: comparar ω_p, γ y τ₀ con las curvas publicadas de Prakapenia–Vereshchagin (2023) y Kluger et al. (texto completo).',
                    'Calcular las condiciones de causalidad de Russo–Townsend con L₁(S,P).',
                    'Buscar sistemáticamente RN/EEH con carga corrida (Uehling) y logarítmico; no se hizo.',
                    'Verificar contra la edición original los datos bibliográficos marcados VERIFY en references.bib (revista/volumen/páginas recordados de memoria) y añadir la afiliación en paper.tex (el autor ya figura: Álvaro López Almeida).',
                    'Incluir los términos de Drummond–Hathrell y restringir 04 a r≳λ_C o rotular como modelo.']
)

cnt = {1: 0, 2: 0, 3: 0}
for r in RESULTS: cnt[r['cat']] += 1
out = dict(meta=dict(date=str(datetime.date.today()), access_limitations='arxiv.org, ar5iv.labs.arxiv.org no resuelven en esta sesión: ningún texto completo leído; solo extractos de WebSearch.',
                     levels=dict(E='extracto de búsqueda esta sesión', M='memoria, no verificado', N='no encontrado en búsqueda superficial', T='texto completo (ninguno)'),
                     counts_by_category=cnt), corrections_to_own_results=[w for w in WEAK if w[0] == 'Error propio detectado'],
           papers=PAPERS, results=RESULTS, weaknesses=[dict(topic=a, detail=b) for a, b in WEAK], publication_hook=HOOK, audit_values=V)
json.dump(out, open('novelty_audit.json', 'w'), indent=1, ensure_ascii=False, default=float)

# ------------------------------- Markdown -------------------------------------------------------------------------
catname = {1: 'Cat.1 Redescubrimiento', 2: 'Cat.2 Refinamiento', 3: 'Cat.3 Síntesis', None: 'Laguna / tensión'}
md = f"""# NOVELTY_AUDIT — Estado del arte y novedad de `DISCOVERY_REPORT.md`
*Generado por `09_make_audit.py` desde `discovery_results.json` (incl. `08_novelty_checks.py`). {out['meta']['date']}.*\n\n> **Actualización tras los cierres teóricos:** P10/P11 (Penrose cargada derivada internamente, con alcance acotado) y 05b (cinética con bloqueo de Pauli). L1 sigue siendo Cat.1; K1 se añade como Cat.2 dependiente de modelo.

## 0. Alcance y límites (léase primero)
- **Ningún texto completo fue leído.** `arxiv.org` y `ar5iv.labs.arxiv.org` no resuelven desde esta sesión (red restringida): hay que permitir esos dominios en *Network access* del entorno para una auditoría de texto completo. Todo lo que sigue se apoya en **extractos de búsqueda** (nivel **E**) o en **memoria** (nivel **M**, no verificada). La *ausencia* de un resultado en una búsqueda superficial (nivel **N**) **no** prueba que no exista.
- Esta auditoría **no** puede dar «rigor absoluto» sobre la novedad: da una clasificación conservadora y una lista de verificaciones pendientes (§7).
- **Corrección de resultados propios** hallada durante la auditoría (ya aplicada en 04, informe y JSON): el margen DEC «−8.33e-4, acotado» era el borde del rango (ε=10⁴⁰); el umbral «1.92» era artefacto de rejilla (exacto **{V['dec_on']:.5f}**). Ahora hay {V['ngates']} gates.

## 1. Veredicto en una tabla
| Resultado | Categoría | En una línea |
|---|---|---|
""" + "".join(f"| **{r['id']}** {r['name']} | {catname[r['cat']]} | {r['equivalent'] if r['equivalent']!='—' else r['ours']} |\n" for r in RESULTS) + f"""
Recuento: Cat.1 = {cnt[1]}, Cat.2 = {cnt[2]}, Cat.3 = {cnt[3]}. **Ningún resultado es, hoy, un descubrimiento físico nuevo defendible**; el valor está en la síntesis cuantitativa, la corrección explícita de la truncación y la verificación formal.

## 2. Comparación artículo por artículo (nivel de verificación entre paréntesis)
| Trabajo | Nivel | Qué se sabe (según lo visto) | Gate(s) del repositorio | Relación |
|---|---|---|---|---|
""" + "".join(f"| {p['ref']} | {p['lvl']} | {p['know']} | {p['gates']} | {catname[p['cat']]} |\n" for p in PAPERS) + """
## 3. Clasificación detallada
"""
for r in RESULTS:
    md += f"""### {r['id']} — {r['name']}  →  **{catname[r['cat']]}**
- **Estructura ya publicada:** {r['evidence']}
- **Ecuación equivalente:** {r['equivalent']}
- **Aporte del repositorio:** {r['ours']}
- **No presentar como:** {r['forbidden']}
- **Cómo citar:** {r['wording']}
- **Salvedades:** {r['caveat']}

"""
md += "## 4. Afirmaciones que NO pueden presentarse como «descubrimiento nuevo» (riesgo de rechazo)\n"
for r in RESULTS:
    if r['cat'] in (1, 2): md += f"- ~~{r['forbidden']}~~ → {r['wording']}\n"
md += "\n## 5. Debilidades que un árbitro señalaría (autocrítica)\n" + "".join(f"{i+1}. **{a}.** {b}\n" for i, (a, b) in enumerate(WEAK))
md += f"""
## 6. El «hook» honesto
**Enunciado.** {HOOK['statement']}

**Por qué es honesto.** {HOOK['why_honest']}

**Contribuciones candidatas (todas condicionadas al §7):**
""" + "".join(f"- {c}\n" for c in HOOK['candidate_claims']) + f"""
**Dónde.** {HOOK['venues']}

**Lo que NO se afirma:** """ + "; ".join(HOOK['not_claimed']) + """.

## 7. Verificaciones pendientes antes de cualquier envío
""" + "".join(f"- [ ] {c}\n" for c in HOOK['before_submission'])
open('NOVELTY_AUDIT.md', 'w').write(md)
print(f"OK: NOVELTY_AUDIT.md ({len(md)} car.) y novelty_audit.json; Cat1={cnt[1]} Cat2={cnt[2]} Cat3={cnt[3]}")
