# Destino del campo eléctrico extremo más allá de la ruptura de la EFT estática
### Descarga de Schwinger, QED no perturbativa y atrapamiento gravitatorio en el sistema Einstein–Maxwell–Euler–Heisenberg

> Informe generado por `07_make_report.py`. **Todo número procede de `discovery_results.json` / `results.json`** (si una clave falta, la generación falla). Cada tesis cita el *gate* que la verifica (`04_*`, `05_*`, `06_*`); el total es 52 gates, todos PASS, de los cuales 10 son controles negativos que **deben** detectar un modelo alterado. Unidades: ħ=c=ε₀=m_e=1, longitudes en λ_C, ε ≡ eE/m_e², G=1/M_Pl².

## Abstract

La línea base (`03_integrate.py`) mostró que la EFT de Euler–Heisenberg truncada falla en ε=1 con ρ ≈ 1.67e-89 M_Pl⁴, muy por debajo de Planck. Aquí se estudia qué ocurre más allá de ese punto con tres herramientas: (A) el Lagrangiano de Euler–Heisenberg a un bucle con **tiempo propio exacto** (parte real e imaginaria), (B) una cinética de pares con retroacción en el régimen de campo localmente constante, y (C) el criterio de atrapamiento de la desigualdad de Penrose cargada. Resultados: (1) la QED a un bucle **no satura** el campo: crece logarítmicamente, `L = (E²/2)[1 − (α/3π) ln ε + …]`, y el polo de Landau está en log₁₀ ε_L = 562, frente a log₁₀ ε_Planck = 44.4: es inaccesible. El perfil interior es Reissner–Nordström con carga corrida (δr₊/r₊ = -1.191e-02 medido vs -1.230e-02 predicho), no el horizonte único tipo Schwarzschild que producen la extrapolación cuártica y Born–Infeld. (2) La descarga neutraliza el campo en **τ_discharge = κ/(e√(eE₀))**, κ = 13.5394, exponente medido -0.50001. (3) En el criterio de Penrose cargada la energía electrostática de un cascarón **se cancela exactamente**: atrapar exige masa no electromagnética, G·M_n ≥ R/2. De ahí se sigue que **la QED no es un censor universal del colapso**: no puede deshacer una superficie atrapada ni crearla desde energía de campo puro; lo que sí hace es destruir la carga neta en τ_dis, siempre más rápido que la caída libre de N masas en reposo para r₀ ≫ r₀* = 9.29e-18 l_P. Se dan el diagrama de fases en (N, r₀) y las leyes de escala. Lo que **no** se resuelve (dinámica tras la inversión del campo, transporte espacial) se declara en §6. **Novedad:** ver `NOVELTY_AUDIT.md`; L1 (masa de Misner–Sharp), la escala de plasma, el crecimiento logarítmico y la violación de DEC (anomalía de traza) son resultados conocidos que aquí se reproducen y verifican formalmente; el valor del trabajo es de síntesis cuantitativa y verificación, no de descubrimiento.

## 1. Línea base y qué pregunta queda abierta

`03_integrate.py` integró el sistema cuártico hacia el centro y se detuvo por ε≥1 (o por escala de variación del campo < λ_C) sin llegar a Planck. Dos hechos se heredan: N_crit = M_Pl²/m_e² = 5.71e+44 (gate **P5**, simbólico: r_Q = r_ε) y ningún horizonte RN cae dentro de la EFT estática si N < N_crit. La pregunta es qué hay en r < r_stop.

*Corrección de rotulado heredado:* `eeh.py:schwinger_rate` es Γ_vac = 2 Im L₁ (decaimiento del vacío, todos los términos de la serie positivos para espinores), **no** el número medio de pares, que es sólo el término n=1 (`qed1loop.pair_rate`). Los criterios de parada del baseline no cambian (usan ε), pero el modelo dinámico usa la tasa correcta. Verificado: gate **N1** (04): |Im L₁(rayo)/serie de Schwinger − 1| < 10⁻⁸.

## 2. Formulación Lagrangiana y cinética

### 2.1 Electrodinámica no lineal estática acoplada (Rama A)
Para **cualquier** L(E) (E=√(g/f)Φ′; X=−2E²), el CAS deriva (`04`, gates **S1–S4**): `T^t_t = T^r_r = −ρ`, `ρ = E L′ − L`, `p_⊥ = T^θ_θ = L`; la primera integral de Maxwell `r² L′(E) = Q/4π`; `y′ = 8πG r² ρ` con f=g=1−y/r (**S2**); y la conservación `ρ′ + (2/r)(ρ+p_⊥) = 0` (**S4**). El control **C1** confirma que un tensor erróneo (p_⊥=L/2) deja residuo −L/r ≠ 0.

Modelos comparados: Maxwell (RN); Euler–Heisenberg cuártico; **QED a un bucle exacto**
`L₁(x) = −(1/8π²) ∫₀^∞ ds/s³ e^{−s}[ x s cot(x s) − 1 + (x s)²/3 ]`, x=eE=ε, evaluado sobre el rayo s=ρe^{−iπ/4} (los polos de cot están en el eje real); y Born–Infeld calibrado `b² = 1/(32c)` (**S5**). Verificaciones: el desarrollo en x reproduce `4cE⁴` con c=α²/90 (**S6**, simbólico) y la correspondencia ε≪1 con el baseline es L₁/(4cE⁴)−1 ∝ ε^2.00 (**N2**); Im L₁ coincide con la serie de Schwinger a ~10⁻¹⁴ (**N1**); la pendiente de campo fuerte es −1/(24π²) con deriva de la constante K=-0.763969 < 10⁻⁹ (**N3**).

### 2.2 Cinética de pares con retroacción (Rama B)
Campo homogéneo, x=eE(t), aproximación de campo localmente constante:

`ẋ = −2e² [ ∫₀ᵗ dt′ Γ(x′)⟨v⟩ + Γ(x)⟨m_⊥⟩/x ]`, `p(t,t′)=∫_{t′}ᵗ x dt″`, `v=p/√(p²+m_⊥²)`, `m_⊥²=1+x s/π` (s~Exp(1), 8 nodos Gauss–Laguerre, **N0**), `Γ(x)=x²/(4π³)e^{−π/x}`.

El primer término es la corriente de conducción y el segundo la de polarización (energía 2m_⊥ de cada par). El CAS verifica: (**D1**) `∇_μT^{μν}_EM = F^ν{}_λ J^λ` en Minkowski esférico con A_t(t,r), A_r(t,r) arbitrarios, de donde `∇_μT^{μν}_matter = −F^ν{}_λ J^λ` si el total se conserva; (**D2**) `d/dt(E²/2 + ρ_pl) = 0` como identidad exacta (Leibniz sobre cohortes). Los controles **C1/C2** (sin j_pol; j_cond con signo invertido) rompen la identidad; **C3/C4** lo detectan en el integrador (|ΔE|/E₀ = 2.9·10⁻² y 1.7·10³³ frente a ~10⁻⁵).

**Validez declarada:** hasta el *primer cero* del campo (las cohortes no se recruzan: sin bloqueo de Pauli ni aniquilación estimulada); LCFA con t_form/τ₀ ≈ e/κ = 0.022 ≪ 1 (**N5**); región ≥ 3 longitudes de formación ⇒ N ≥ 9/α = 1233.

### 2.3 Criterio de atrapamiento (Rama C)
Para una superficie atrapada exterior de radio R con toda la carga en su interior, la desigualdad de Penrose cargada exige `G M_ADM ≥ (R + r_Q²/R)/2`, r_Q² = G Q²/4π (Malec–Ó Murchadha; Khuri–Weinstein–Yamada; *citada, no demostrada aquí*). El CAS verifica que RN la **satura** (**P1**) y que para un cascarón `G U_C(R) = r_Q²/(2R)` (**P2**), de modo que `G(M_n+U_C) − (R+r_Q²/R)/2 = G M_n − R/2` (**P3**). Con N electrones confinados (masa en reposo + energía de Fermi ultrarrelativista) el radio de atrapamiento es la raíz positiva `r_s(N)` de `r = 2G(N + (3/4)(9π/4)^{1/3}N^{4/3}/r)` (**P4**).

## 3. Leyes de escala y resultados deducidos

**L1 (Cancelación de Coulomb).** El atrapamiento de un cascarón cargado no depende de Q: exige G M_n ≥ R/2 (P2, P3). *Corolario:* una superficie atrapada formada sólo con energía de campo es imposible (`verdict.pure_field_trapped_possible = False`). La razón carga/masa de N electrones es z² = α M_Pl² = 4.17e+42 veces la extremal, y una carga unitaria extremal pesaría 2.04e+21 m_e.

**L2 (Tiempo de descarga).** `τ_discharge = κ /(e√(eE₀))` con κ∞ = 13.5394 (dispersión 2.6e-05 en ε₀∈[10⁶,10¹⁰]); exponente -0.50001 (**N4**). Convergencia: orden observado 1.98 con h, h/2, h/4 (**N1**); método independiente (trapecio implícito) τ₀ = 1.416579 frente a 1.416579 (**N2**); la energía discreta se conserva a 1.3e-05 (suelo ~10⁻⁵ en el cruce; **N3**). Con E₀ = Q/(4πr₀²) en la superficie, ε₀ = Nα/r₀² y `τ_discharge = κ r₀/(α√(4πN))`: la carrera con la luz es **independiente de r₀**; τ_dis = r₀ en N* = κ²/(4πα²) = 2.74e+05 (**P7**).

**L3 (Carrera con el colapso).** Para el colapso libre de N masas en reposo, τ_ff = (π/2)√(r₀³/2GM) (**P8**, integral simbólica), y
`τ_discharge / τ_ff = √(r₀*/r₀)`, `r₀* = 2κ²/(π³α²M_Pl²) = 3.89e-40 λ_C = 9.29e-18 l_P`
(**P9**; **la dependencia en N se cancela exactamente**). Como r₀* ≪ l_P, la descarga es más rápida que la caída libre para todo r₀ físico. En la frontera de atrapamiento r₀=r_s con energía de Fermi incluida, la descarga gana (τ_dis<τ_ff) sólo para N ≳ 1.11e+05; por debajo, el colapso es más rápido, **sin que esto cambie el veredicto** (L1: la carga no entra en el criterio).

**L4 (Rama A: campo fuerte no perturbativo).** `L = (E²/2)[1 − (α/3π)(ln ε + 3K)]` (**S7**, simbólico), polo en `ln ε_L = 3π/α − 3K` (**S8**): log₁₀ε_L = 561.9 vs log₁₀ε_Pl = 44.4, es decir ln ε_L es 12.7 veces ln ε_Pl: **el polo de Landau es inaccesible a un bucle** (ver el límite de Ritus–Narozhny en §6). WEC y NEC se cumplen hasta ε=10⁴⁰ (min ρ/E² = 0.465, min ED/E² = 0.930, **N7**); la SEC (L>0) vale hasta el polo (**N8**); la **DEC se viola** para ε > 1.8822 (umbral exacto = máximo de L₁/x², identidad simbólica **S10**). El margen relativo **no está acotado**: tiende a −α/(3π) = -7.743e-04 y crece despacio con ln ε (-7.875e-04 en ε=10¹⁰, -8.172e-04 en 10³⁰) hasta divergir en el polo (**N9**). Es la **anomalía de traza** a un bucle (coeficiente β de la QED), física conocida, no un efecto nuevo; el valor -8.328e-04 de versiones previas era el borde de mi rango (ε=10⁴⁰). D′(E) > 0 siempre: no hay pérdida de hiperbolicidad (**N4**).

## 4. Resultados numéricos

### 4.1 Perfil interior del campo y de la métrica (incógnita a)
`L/(E²/2)`: en ε=10 cuártico 1.0052 vs QED 0.99967; en ε=10⁶, 5.162e+07 vs 0.99108; en ε=10³⁰ la QED sigue en 0.9483. La truncación α² **sobreestima catastróficamente** la no linealidad.

Métrica con masa M=1.5·M_ext y N=10⁶ (r_Q=3.575e-18):

| Modelo | horizontes (r₋, r₊) |
|---|---|
| RN (analítico) | 1.36569e-18, 9.36055e-18 |
| QED a un bucle | 1.48259e-18, 9.24904e-18 |
| EH cuártico extrapolado | único: 1.07263e-17 |
| Born–Infeld calibrado | único: 1.07263e-17 |

El valor único del cuártico y de BI es 2GM (Schwarzschild): ambos eliminan la energía de campo interior. La QED reproduce RN con carga corrida: λ=E/E_Coulomb(r₊)=1.0710 y δr₊/r₊ medido -1.191e-02 vs predicho -1.230e-02 (**N6**). Con N=10²⁰ el patrón se repite (`branch_A.profiles.1e+20`). Resultado: **la no linealidad QED estática no censura el horizonte** (lo desplaza ~1%); la respuesta no puede venir de la estática.

### 4.2 Escalas por número de cargas
| N | r_ε=√(Nα) | r_Q (RN) | r_s (atrapamiento) | ε(r_s) | τ_dis/r0 |
|---|---|---|---|---|---|
| 1 | 0.0854 | 3.58e-24 | 7.1e-23 | 1.4e+42 | 523 |
| 1e+03 | 2.7 | 3.58e-21 | 7.1e-21 | 1.4e+41 | 16.6 |
| 1e+06 | 85.4 | 3.58e-18 | 7.1e-19 | 1.4e+40 | 0.523 |
| 1e+12 | 8.54e+04 | 3.58e-12 | 7.1e-15 | 1.4e+38 | 0.000523 |
| 1e+20 | 8.54e+08 | 0.000358 | 1.53e-09 | 3.1e+35 | 5.23e-08 |
| 1e+30 | 8.54e+13 | 3.58e+06 | 0.0071 | 1.4e+32 | 5.23e-13 |
| 5.71e+44 | 2.04e+21 | 2.04e+21 | 4.89e+07 | 1.7e+27 | 2.19e-20 |
| 1e+48 | 8.54e+22 | 3.58e+24 | 7.1e+09 | 1.4e+26 | 5.23e-22 |
| 1e+60 | 8.54e+28 | 3.58e+36 | 7.12e+17 | 1.4e+22 | 5.23e-28 |

(`06`; r_s con N electrones degenerados; τ_dis/r₀ sólo es aplicable con N ≥ 1233.)

### 4.3 Diagrama de fases (incógnita c)
![Diagrama de fases](phase_diagram.png)

Fases en el plano (N, r₀) (fracción del plano explorado entre paréntesis): (0) EFT estática, ε<1 (0.310); **(i) descarga de Schwinger** en espacio casi plano, r_P < r₀ < r_ε (0.281); **(ii) atrapamiento/horizonte**, r₀ ≤ r_s(N) (0.087); **(iii) trans-Planck**, ρ_Coul ≥ M_Pl⁴ o r₀<10 l_P (0.322). Umbrales: N_P=(8πα)^{−1/2}=2.335 (**P6**: ρ(r_Q)/M_Pl⁴ = 1/(8παN²)); atrapamiento sub-planckiano desde N≈15 (**N3**); N_crit=5.71e+44; un BH clásico entra en la EFT estática (r_s>r_ε) sólo para N ≳ α M_Pl⁴/4 = 5.94e+86 (**N5**). La fase (ii) está **siempre por debajo** de r_ε (ε≫1) y por encima de Planck: el atrapamiento ocurre con densidad sub-planckiana y campo eléctrico supercrítico.

## 5. Rutas descartadas (con su evidencia)
1. **Extrapolar el EH cuártico a ε≫1.** E∼r^{−2/3}, ρ∼r^{−8/3}, masa interior finita (**S9**): artefacto; la QED a un bucle da ρ∼r^{−4}.
2. **Born–Infeld como completación de la QED.** Calibrado a b²=1/32c satura en ε=69.6, pero ya a ε=10 el signo de la desviación respecto a Maxwell es opuesto al de la QED (BI 1.0052 > 1; QED 0.99967 < 1). No es la completación de la QED; se conserva sólo como contraste.
3. **Censura por saturación estática.** Refutada por 4.1 (desplazamiento ~1%).
4. **Densidad de energía de campo puro como fuente de horizonte.** Imposible por L1 (cancelación exacta).

## 6. Límites UV, alcance y problemas abiertos (sin extrapolaciones silenciosas)
- **Frontera de predicción:** el polo de Landau y la densidad de Planck no se cruzan en ningún cálculo; el criterio ρ ≥ M_Pl⁴ / r<10 l_P define la fase (iii) y ahí no se afirma nada.
- **No resuelto:** (a) la dinámica **tras el primer cero** del campo (oscilaciones de plasma, bloqueo de Pauli con ocupación f→1 y aniquilación estimulada): requiere la ecuación cinética cuántica completa con renormalización; (b) transporte espacial (esfera con perfil, hidrodinámica o Vlasov 1+1) y choques/cruce de flujos; (c) el modelo es plano y homogéneo: la métrica dependiente del tiempo sólo entra por el criterio de atrapamiento y la contabilidad de energía; (d) la hipótesis «la descarga neutraliza completamente la carga neta por debajo de r_neutral ≈ 0.085 λ_C» es una **estimación no verificada**; (e) κ tiene incertidumbre sistemática O(1) por LCFA y por usar la tasa de pares de campo constante; los *exponentes* son robustos, los prefactores no.
- **Ritus–Narozhny (omisión corregida tras la auditoría):** la conjetura de que cada lazo aporta αχ^{2/3} (χ≳α^{-3/2}≈1.6·10³ en campo cruzado constante) implica que todo resultado «a un bucle» para ε≫10³ —incluidos WEC/NEC/SEC hasta 10⁴⁰, el polo de Landau y el barrido de 05 hasta ε=10¹⁰— puede alterarse a orden superior; su aplicabilidad a campo eléctrico puro no está verificada.
- **Expansión derivativa:** los perfiles de 04 en r∼r_Q≪λ_C (ℓ≈r_Q/2≪λ_C) violan el criterio derivativo del baseline y omiten los términos de curvatura de Drummond–Hathrell (α(λ_C/L)²); son soluciones *formales* del modelo de campo localmente constante, no predicciones de la QED.
- **Hipótesis citadas, no demostradas en el repositorio:** desigualdad de Penrose cargada (requiere toda la carga dentro de la superficie atrapada) y persistencia de superficies atrapadas bajo NEC. Que un agujero negro cargado pierde carga por Schwinger es el resultado de Gibbons (1975); aquí sólo se comprobó la consistencia cualitativa (τ_dis ∼ r₀ para N≈10⁵–10⁶).
- La fase (ii) usa un modelo de materia (electrones degenerados) que es un *parámetro*: otro M_n desplaza r_s, no la estructura lógica (columna `r_s_ext` en `branch_C.table`).

## 7. Conclusión
Más allá de r_stop no aparece ninguna singularidad ni saturación de origen QED en el régimen calculable: la no linealidad estática crece sólo como α ln ε y el polo de Landau queda 12.7 veces más lejos (en ln ε) que Planck. La pregunta «¿es la QED un censor cuántico universal del colapso?» tiene respuesta **negativa**, con matices precisos: (1) una superficie atrapada exige masa no electromagnética G M_n ≥ R/2; la energía electrostática no cuenta (cancelación exacta, L1); (2) la descarga de Schwinger convierte el campo en un plasma de pares en τ_dis = κ/(e√(eE₀)), conservando energía-impulso, y por tanto no puede reducir M_n ni deshacer una superficie atrapada; (3) lo que **sí** hace la QED es eliminar la carga neta que, por ser z² ≈ 4·10⁴² veces superextremal, impediría por repulsión el colapso de la materia cargada: actúa como *facilitador*, más rápido que la caída libre (L3). La condición real de formación de un horizonte a escalas sub-Compton es de masa, no de carga: N m_e ≳ r₀/2G, es decir N ≳ M_Pl²/2 para r₀=λ_C, que es el N_crit de la línea base reinterpretado.

## Apéndice A. Trazabilidad (gate → script)
| Gate | Script | Estado |
|---|---|---|
| S1 T^mu_nu generico L(E) | 04_*.py | PASS |
| S2 Einstein tt -> y' | 04_*.py | PASS |
| S3 Maxwell | 04_*.py | PASS |
| S4 conservacion generica | 04_*.py | PASS |
| C1 control | 04_*.py | PASS |
| S5 BI -> cuartico (b^2=1/32c) | 04_*.py | PASS |
| S6 tiempo propio -> cuartico | 04_*.py | PASS |
| S7 coeficiente del running | 04_*.py | PASS |
| S8 polo de Landau (simbolico) | 04_*.py | PASS |
| N1 Im L1 = serie de Schwinger | 04_*.py | PASS |
| N2 correspondencia ε≪1 -> EH cuartico (baseline) | 04_*.py | PASS |
| N3 pendiente logaritmica | 04_*.py | PASS |
| N4 hiperbolicidad D'(E)>0 (QED 1-bucle) | 04_*.py | PASS |
| N5 RN recuperado (Maxwell) | 04_*.py | PASS |
| N6 horizonte QED = RN con carga corrida | 04_*.py | PASS |
| S9 exponente cuartico E~r^-2/3 | 04_*.py | PASS |
| S10 identidad DEC | 04_*.py | PASS |
| N9 DEC | 04_*.py | PASS |
| N7 WEC/NEC cumplidas hasta Planck | 04_*.py | PASS |
| N8 SEC (L>0) válida hasta Landau | 04_*.py | PASS |
| C2 control | 04_*.py | PASS |
| C3 control | 04_*.py | PASS |
| D1 ∇_μT^{μν}_EM = F^{ν}{}_λ J^λ (esferico, A_t(t,r),A_r(t,r)) | 05_*.py | PASS |
| D2 conservacion de energia campo+plasma | 05_*.py | PASS |
| C1 control | 05_*.py | PASS |
| C2 control | 05_*.py | PASS |
| D3 Gamma ~ x^2/(4pi^3) (ε→∞) | 05_*.py | PASS |
| N0 cuadratura Gauss-Laguerre | 05_*.py | PASS |
| N1 convergencia h,h/2,h/4 (Heun) | 05_*.py | PASS |
| N2 metodo independiente (trapezoide implicito) | 05_*.py | PASS |
| N3 conservacion de energia (discreta) | 05_*.py | PASS |
| N4 ley de escala τ0 ∝ ε0^(-1/2) | 05_*.py | PASS |
| N5 LCFA válida | 05_*.py | PASS |
| C3 control numerico | 05_*.py | PASS |
| C4 control numerico | 05_*.py | PASS |
| P1 RN satura Penrose cargada | 06_*.py | PASS |
| P2 cancelacion Coulomb (cascaron) | 06_*.py | PASS |
| P3 criterio reducido | 06_*.py | PASS |
| P4 radio de atrapamiento r_s(N) | 06_*.py | PASS |
| P5 N_crit = M_Pl^2 | 06_*.py | PASS |
| P6 densidad en r_Q | 06_*.py | PASS |
| P7 N* (τ_dis = r0) | 06_*.py | PASS |
| P8 caida libre OS | 06_*.py | PASS |
| C1 control | 06_*.py | PASS |
| P9 ley de la carrera | 06_*.py | PASS |
| N1 N_crit numerico | 06_*.py | PASS |
| N2 r_s con N electrones | 06_*.py | PASS |
| N3 BH por encima de la region trans-planckiana | 06_*.py | PASS |
| N4 N_P | 06_*.py | PASS |
| N5 BH clasico dentro de la EFT estatica | 06_*.py | PASS |
| C2 control | 06_*.py | PASS |
| C3 control | 06_*.py | PASS |

## Apéndice B. Referencias (leídas como extractos de búsqueda; no verificadas en texto completo)
- Cinética cuántica del mecanismo de Schwinger y límite de Markov: [Schmidt et al., hep-ph/9809227](https://arxiv.org/pdf/hep-ph/9809227); [hep-ph/9712377](https://arxiv.org/pdf/hep-ph/9712377); [Pair production in a strong electric field, hep-ph/0311293](https://arxiv.org/abs/hep-ph/0311293). Bloqueo de Pauli en campos sobrecríticos: [arXiv:2306.12956](https://arxiv.org/pdf/2306.12956).
- Pérdida espontánea de carga de agujeros negros: Gibbons, Commun. Math. Phys. 44 (1975) 245; [Gabriel, gr-qc/0010103](https://arxiv.org/pdf/gr-qc/0010103).
- Desigualdad de Penrose cargada y superficies atrapadas en simetría esférica: [Malec–Ó Murchadha, gr-qc/9401014](https://arxiv.org/pdf/gr-qc/9401014); [Khuri–Weinstein–Yamada](https://www.math.stonybrook.edu/~khuri/Khuri_Weinstein_Yamada_Cont_Math.pdf); [Kunduri, arXiv:2409.07639](https://arxiv.org/html/2409.07639v1).
- Heisenberg–Euler, Schwinger (1951), Born–Infeld, Oppenheimer–Snyder: citados de memoria (no consultados en esta sesión).
