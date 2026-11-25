# Roteiro de Desenvolvimento — Otimização de Trajetórias de Múltiplos Sobrevoos em Asteroides (AG Real & Algoritmo Memético)

Este documento estabelece o roteiro técnico completo, incremental e verificável para o desenvolvimento em **C++17/20** (CMake, compilação com `-O3`) de um **Algoritmo Genético de Codificação Real** e sua extensão **Memética (hibridizada com Busca Local Incremental)** aplicados à **Etapa 2** (otimização de datas contínuas sobre sequências fixas) dos benchmarks **GTOC4 Impulsivo** e **GTOC11**, baseados em Zhang et al. (2025).

> [!IMPORTANT]
> **Diretórios Somente Leitura vs. Implementação:**
> - `Fontes/` e `Multi_Flyby_Dynamic_Programming/` são **estritamente somente leitura**.
> - Todo código-fonte, `CMakeLists.txt`, build (`Codigo/build/`), testes, scripts e saídas experimentais devem residir exclusivamente em `Codigo/` *(nota: existe atualmente no disco uma pasta vazia acentuada `Código/`; na Fase 0 ela deve ser padronizada para `Codigo/`, sem acento, para compatibilidade com CMake e toolchains Linux)*.

---

## 1. Mapeamento Real do Workspace, Fontes, Constantes e Benchmarks

### 1.1 Arquivos Inspecionados em `Fontes/` (somente leitura)

| Arquivo | Descrição e Papel no Projeto |
|---|---|
| `Fontes/2508.02904v1.pdf` | Zhang, Guo, Wu, Baoyin, Li & Topputo (2025), *"Global Optimality in Multi-Flyby Asteroid Trajectory Optimization: Theory and Application Techniques"* (36 páginas). Referências-chave: **Seção II** (formulação $P_0/P_1/P_2$, Eqs. 1–12), **Seção IV.E.2** (*Bi-Impulse*, Alg. 1, Eq. 34), **Seção V.A** (GTOC4 Impulsivo, Tabelas 1 e 2), **Seção V.B** (GTOC11, Eq. 35, Tabelas 3 e 4) e **Seção V.C** (parâmetros do GTOC4, Tabela 5). |
| `Fontes/Enunciado do trabalho 2 DCC136 - 2026-3.pdf` | Enunciado do Trabalho Final de DCC136 (*Inteligência Computacional — UFJF 2026.3*, entrega final em 07/12/2026). Exige: modelagem computacional própria, implementação da metaheurística pelo grupo (sem delegar a bibliotecas prontas), definição e justificativa de parâmetros, múltiplas execuções com diferentes sementes aleatórias, análise crítica frente à literatura, código limpo/documentado e **proibição de incluir código-fonte no relatório ou apresentações**. |

### 1.2 Arquivos Inspecionados em `Multi_Flyby_Dynamic_Programming/` (somente leitura)

| Caminho Relativo | Conteúdo Real e Observações Críticas para Reaproveitamento |
|---|---|
| `code/include/Constant.h` | Constantes matemáticas (`DPI`, `D2PI`, `D2R`, `R2D`, `EPSILON = 1.0e-14`) e físicas (`JD2Second = 86400.0`, `g0 = 9.80665`, `muSun = 1.32712440018e20`, `AU = 1.49597870700e11`). **Atenção:** `GTOC4_Problem.h` e `GTOC11_problem.h` definem `AU = 1.49597870691e11` m; deve-se usar exatamente `1.49597870691e11` m para reprodução numérica fiel dos benchmarks. |
| `code/include/OrbitMath.h` | Utilitários inline de álgebra vetorial/matricial (`V_Add`, `V_Minus`, `V_Dot`, `V_Cross`, `V_Norm2`, `V_Vers`, `M_Multi`, `M_Tranpose`, `Sign`, `NiceAngle`). |
| `code/include/OrbitFun.h`<br>`code/include/OrbitFun.cpp` | Conversão entre anomalias (`E2f`, `E2M`, `f2E`, `M2E` via esquema iterativo de 3ª ordem com $\varepsilon = 10^{-14}$, `MaxIter = 100`), propagação kepleriana analítica (`f0dt2ft`, `coe02coef`, `rv02rvf`) e conversões (`coe2rv`, `rv2coe`). Nota: `f0ft2dt` está declarada no `.h`, mas não implementada no `.cpp`. |
| `code/include/Lambert.h`<br>`code/include/Lambert.cpp` | Solucionador de Izzo (2004) traduzido por Jiang FH (`x2tof`, `lambert`), busca multi-revolução (`LambEval`, varrendo $N \in [0, N_{\max}]$ com $N_{\max} = \lfloor t / (2\pi\sqrt{a_m^3/\mu}) \rfloor$ e ramos `branch = 0, 1`) e classe `lambert_solver` (`solve_single`, `solve_multi`, `get_dv1`, `get_dv2`, `get_Mdv1`, `get_Mdv2`). |
| `code/input/GTOC4_data.txt` | **1.437 linhas**: 2 linhas de cabeçalho iniciadas com `#` + **1.436 asteroides** (índices 0-based `0..1435`, última linha sem `\n` final). Colunas separadas por espaço: `Name`, `Epoch (MJD)`, `a (AU)`, `e`, `i (deg)`, `LAN (deg)`, `arg. periap. (deg)`, `mean anomaly (deg)`. **Fato crítico:** 1.435 asteroides têm época `54800.0` MJD, mas **1 asteroide tem época `54745.918153` MJD**. O parser deve ler a época individual da linha e propagar para `54800.0` MJD (ou armazenar `epoch_mjd` por corpo). |
| `code/input/Candidate_Asteroids.txt` | **83.453 linhas**, **sem linha de cabeçalho**, separadas por TAB (`\t`). 9 colunas: `id (1..83453)`, `epoch (MJD, todos = 59396)`, `a (AU)`, `e`, `i (deg)`, `Omega (deg)`, `omega (deg)`, `M (deg)`, `mass (kg)`. |
| `code/src/GTOC4_Problem.h`<br>`code/src/GTOC4_Problem.cpp` | Constantes GTOC4, `read_asteroid_data_GTOC4`, `getcoe`/`getrv` (Terra = `-3`, asteroides = `0..1435`), `compute_v_nominal`, `chimpionship_result` (sequência impulsiva de **50 corpos** convertida de 1-based para 0-based com Terra = `-3` e épocas `T_GTOC4` ordenadas por `std::sort`) e `jena_result` (sequência low-thrust de 51 corpos — não utilizada na Etapa 2 impulsiva). |
| `code/src/GTOC11_problem.h`<br>`code/src/GTOC11_problem.cpp` | Constantes SI e normalizadas do GTOC11, `read_data_GTOC11`, `getcoe_GTOC11`/`getrv_GTOC11` (Terra = `-2`, asteroides = `0..83452` indexados diretamente em `data_GTOC11`), `GTOC11_sequence_process` (`rawData`, `Sequence_GTOC11_multi`, `T_GTOC11_multi`, `Maximum_T_GTOC11` para as **10 missões**), e método cônico `rendezvous2flyby` (`calc_transF2X`, `verifyFinC`, `find_minimum`, `calc_eqnX`). |
| `code/src/DP_flyby_2impulse.h/.cpp`<br>`code/src/DP_flyby_2impulse_user_function.h/.cpp`<br>`code/src/main.cpp` | Referência de cálculo de $\Delta V$ por estágio (`compute_rv_dv`), definição de janelas/ToF (`possible_t_values`, `check_t`) e formatação de saída (`print_optimal_path`). |
| `result_solution/DP_result_for_orignal_GTOC4.txt` | **Esclarecimento crítico:** contém apenas a solução de **baixo empuxo** do GTOC4 original (`jena_result`, 51 estágios, massa final `519,9 kg`), **não** a solução impulsiva. |
| `code/bin/DP_result.txt` | Contém as trajetórias **impulsivas** completas geradas pelo código de referência:<br>- Linhas 1–38: GTOC11 Missão 1 (`Total DV = 13944.823 m/s`).<br>- Linhas 39–76: execução alternativa da Missão 1 (`14261.473 m/s`).<br>- Linhas 77–456: GTOC11 Missões 2 a 10.<br>- Linhas 457–664: GTOC4 Impulsivo com épocas do campeão (`Total DV = 26167.447 m/s`, `Mass = 616.325 kg`). |

### 1.3 Constantes Físicas Exatas e Efemérides da Terra

#### Caso 1 — GTOC4 Impulsivo
- $\mu_{\text{Sol}} = 1.32712440018 \times 10^{20}\;\text{m}^3/\text{s}^2$ (`Mu_GTOC4`)
- $\text{AU} = 1.49597870691 \times 10^{11}\;\text{m}$ (`AU_GTOC4 * 1000.0`)
- $g_0 = 9.80665\;\text{m/s}^2$, $I_{\text{sp}} = 3000\;\text{s}$, $m_0 = 1500.0\;\text{kg}$, $m_{\min} = 500.0\;\text{kg}$
- $\text{Dia} = 86400.0\;\text{s}$, $\text{Ano} = 365.25\;\text{dias}$, duração máxima $t_{M-1} - t_0 \le 10 \times 365.25 = 3652.5\;\text{dias}$
- $V_{\infty,\text{Terra},\max} = 4000.0\;\text{m/s}$
- `Initial_MJD_GTOC4 = 57023.0` (referência de conversão $t_{\text{seg}} = (\text{MJD} - 57023.0) \cdot 86400.0$)
- **Terra (`ID = -3`) na época `54000.0 MJD`:**
  - $a = 0.999988049532578 \times \text{AU}$ (m)
  - $e = 1.671681163160 \times 10^{-2}$
  - $i = 0.8854353079654 \times 10^{-3}\;\text{deg}$
  - $\Omega = 175.40647696473\;\text{deg}$
  - $\omega = 287.61577546182\;\text{deg}$
  - **Anomalia Verdadeira** $f = 4.4635844709062136\;\text{rad}$ *(já convertida de $M$ no código original — **não** aplicar `M2E`/`E2f` novamente)*
- **Asteroides do GTOC4 (`ID = 0..1435`):** propagados a partir de sua época individual (`54800.0` ou `54745.918153` MJD) lida do arquivo.

#### Caso 2 — GTOC11
- $\mu_{\text{Sol}} = 1.32712440018 \times 10^{20}\;\text{m}^3/\text{s}^2$ (`GTOC11_mu`)
- $\text{AU} = 1.49597870691 \times 10^{11}\;\text{m}$ (`GTOC11_AU`)
- Unidades normalizadas:
  - $L_{\text{unit}} = \text{AU}$
  - $T_{\text{unit}} = \sqrt{L_{\text{unit}}^3 / \mu_{\text{Sol}}} \approx 5022642.890913\;\text{s} \approx 58.132441\;\text{dias}$
  - $V_{\text{unit}} = L_{\text{unit}} / T_{\text{unit}} \approx 29784.691833\;\text{m/s}$
- `GTOC11_MJD_Init = 95739.0`, `GTOC11_MJD_End = 103044.0` (referência $t_{\text{seg}} = (\text{MJD} - 95739.0) \cdot 86400.0$)
- $V_{\infty,\text{Terra},\max} = 6000.0\;\text{m/s}$, velocidade relativa máxima de sobrevoo $v_{\text{rel,flyby},\max} = 2000.0\;\text{m/s}$ (`r = (2.0 - 1.0e-12) * 1000.0` m/s em `rendezvous2flyby`)
- **Terra (`ID = -2`) na época `59396.0 MJD`:**
  - $a = 9.998012770769207 \times 10^{-1}\;\text{AU}$
  - $e = 1.693309475505424 \times 10^{-2}$
  - $i = 3.049485258137714 \times 10^{-3}\;\text{deg}$
  - $\Omega = 166.2869706216879\;\text{deg}$
  - $\omega = 297.8214889887391\;\text{deg}$
  - **Anomalia Média** $M = 175.7352290983351\;\text{deg}$ *(converter para rad e depois via `M2E` $\to$ `E2f`)*
- **Asteroides do GTOC11 (`ID = 0..83452`):** nota crítica de indexação! Embora `Candidate_Asteroids.txt` tenha coluna 1 com valores `1..83453`, `GTOC11_problem.cpp` lê sequencialmente em `data_GTOC11[0..83452]` e as sequências em `rawData` usam diretamente índices 0-based em `data_GTOC11` (ex.: `74927`, `83419`).

### 1.4 Sequências Fixas, Janelas de Busca `[L_k, U_k]` e Valores de Referência (BKS)

#### Sequência e BKS do GTOC4 Impulsivo
- **Comprimento:** $M = 50$ corpos (`-3` Terra + 48 sobrevoos + 1 encontro final no asteroide `366` = `2000SZ162`). Sequência 0-based completa extraída de `chimpionship_result`:
  `{-3, 1057, 1125, 1328, 1181, 937, 1077, 1289, 1000, 601, 1376, 1253, 954, 1090, 1022, 752, 151, 958, 1337, 1365, 800, 1356, 1245, 1329, 666, 926, 735, 1083, 234, 815, 987, 132, 833, 222, 939, 1293, 709, 1238, 1379, 248, 1118, 475, 651, 970, 743, 1100, 555, 1385, 1173, 366}`.
- **Janelas de busca `[L_k, U_k]` (em MJD, sem vazamento da solução BKS):**
  - Seguindo exatamente `DP_flyby_2impulse_user_function.cpp`: $T_{\text{win},0} = 58676.4$, $T_{\text{win},49} = 62261.25$, $t_{\text{mid}} = 60468.825\;\text{MJD}$, envelope global $[T_{\min}, T_{\max}] = [58642.575,\; 62295.075]\;\text{MJD}$ (10 anos = $3652.5$ dias).
  - Centro linear uniforme: $t_{\text{center},k} = T_{\min} + (k + 0.5)\frac{T_{\max} - T_{\min}}{49}$.
  - Janela por estágio: $L_k = \max(T_{\min},\; t_{\text{center},k} - 480.0)$, $U_k = \min(T_{\max},\; t_{\text{center},k} + 480.0)$ (correspondente a $\pm 15 \times 32 = 480$ dias).
  - **Verificação:** tanto o campeão (desvio máximo de $138{,}28\;\text{d}$) quanto a BKS da Tabela 2 (desvio máximo de $136{,}13\;\text{d}$) pertencem estritamente a $[L_k, U_k]$.
  - Limites de tempo de voo por arco: $\text{ToF}_{\min} = 1.0\;\text{dia}$, $\text{ToF}_{\max} = 365.25\;\text{dias}$.
- **Valores de Referência (GTOC4 Impulsivo):**
  - Campeão com datas `T_GTOC4` (`DP_result.txt`): $\Delta V = 26.167{,}447\;\text{m/s}$, $m_f = 616{,}325\;\text{kg}$.
  - **BKS Artigo (Tabelas 1 e 2):** $\Delta V = 25.210{,}4\;\text{m/s}$, $m_f = 636{,}7\;\text{kg}$.

#### Sequências, Janelas e BKS do GTOC11 (10 Missões)
- **Sequências e $T_{\max,k}$:** extraídos literalmente de `rawData` em `GTOC11_problem.cpp` (onde `rawData[m][0]` = $M$, `rawData[m][1..M]` = IDs iniciando em `-2`, e `rawData[m][M+1..2M]` = `Maximum_T` em unidades normalizadas de tempo $T_{\text{unit}}$).
- **Janelas de busca `[L_k, U_k]` (em MJD):**
  - Limite global inferior: $T_{\min} = 95739.0\;\text{MJD}$.
  - Limite superior individual de Fase 2: $T_{\max,k} = \text{Maximum\_T}[m][k] \cdot T_{\text{unit}} / 86400.0\;\text{MJD}$.
  - Centro linear uniforme: $t_{\text{center},k} = T_{\min} + (k + 0.5)\frac{T_{\max,M-1} - T_{\min}}{M - 1}$.
  - **Atenção (descoberta verificada no código de referência):** enquanto a grade inicial fixa do DP usa $\pm 480\;\text{d}$ e desloca o tubo adaptativo entre iterações, as soluções ótimas das Missões 1 e 3 chegam a desviar $501{,}13\;\text{d}$ e $513{,}08\;\text{d}$ de $t_{\text{center},k}$. Para evitar cortar regiões viáveis ótimas no AG contínuo sem usar as datas da solução de referência, adotar meia-largura de **$\pm 600{,}0\;\text{dias}$**:
    $$L_k = \max(T_{\min},\; t_{\text{center},k} - 600.0), \qquad U_k = \min(T_{\max,k},\; t_{\text{center},k} + 600.0)$$
  - Limites de tempo de voo por arco: $\text{ToF}_{\min} = 1.0\;\text{dia}$, $\text{ToF}_{\max} = 730.5\;\text{dias}$ (2 anos).

| Missão | Corpos ($M$) | Alvos ($M-1$) | $T_{\max,0}$ (MJD) | $T_{\max,M-1}$ (MJD) | Winner Bi-Impulse (m/s) | `DP_result.txt` (m/s) | **BKS Artigo Tab. 3 (m/s)** |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 36 | 35 | 97763,00 | 101416,03 | 14.325,9 | 13.944,823 | **13.944,8** |
| 2 | 43 | 42 | 97415,00 | 101890,49 | 17.199,7 | 17.147,721 | **17.131,1** |
| 3 | 40 | 39 | 97423,00 | 101644,02 | 17.330,1 | 16.764,969 | **16.704,1** |
| 4 | 42 | 41 | 97515,99 | 101665,00 | 18.508,0 | 18.381,558 | **18.351,0** |
| 5 | 41 | 40 | 97335,49 | 101729,48 | 16.680,7 | 16.659,698 | **16.033,9** |
| 6 | 40 | 39 | 97135,00 | 101704,52 | 19.824,3 | 19.728,864 | **19.721,8** |
| 7 | 38 | 37 | 97093,51 | 101636,44 | 17.606,0 | 17.376,584 | **17.352,9** |
| 8 | 41 | 40 | 97521,50 | 101914,52 | 16.342,2 | 16.335,086 | **16.331,0** |
| 9 | 38 | 37 | 97720,01 | 101783,50 | 15.158,4 | 15.061,754 | **15.025,1** |
| 10 | 39 | 38 | 97203,50 | 101791,97 | 16.168,8 | 16.142,227 | **16.126,6** |

---

## 2. Modelagem Física, Arquitetura do Algoritmo e Estrutura de `Codigo/`

### 2.1 Função Objetivo e Tratamento de Restrições

Seja $X = [t_0, t_1, \dots, t_{M-1}]$ o vetor contínuo de datas em **MJD** sobre uma sequência fixa $[b_0, b_1, \dots, b_{M-1}]$, com $M-1$ arcos de transferência kepleriana resolvidos por `solve_lambert_arc` (multi-revolução, menor $\Delta v$ de transferência). Para cada arco $a \in \{0, \dots, M-2\}$ entre $(b_a, t_a)$ e $(b_{a+1}, t_{a+1})$, obtêm-se $\mathbf{v}_{\text{dep}}[a]$ (velocidade heliocêntrica da nave na saída de $b_a$) e $\mathbf{v}_{\text{arr}}[a]$ (velocidade heliocêntrica da nave na chegada a $b_{a+1}$).

#### Caso 1: GTOC4 Impulsivo ($M = 50$)
1. **Nó de Partida ($k = 0$, Terra):**
   $$v_{\infty,0} = \|\mathbf{v}_{\text{dep}}[0] - \mathbf{V}_{\text{Terra}}(t_0)\|_2, \qquad \Delta v_0 = \max(0,\; v_{\infty,0} - 4000.0)$$
2. **Nós de Sobrevoo Intermediários ($k = 1 \dots M-2$):**
   $$\Delta v_k = \|\mathbf{v}_{\text{dep}}[k] - \mathbf{v}_{\text{arr}}[k-1]\|_2$$
3. **Nó de Encontro Final ($k = M-1$):**
   $$\Delta v_{M-1} = \|\mathbf{V}_{b_{M-1}}(t_{M-1}) - \mathbf{v}_{\text{arr}}[M-2]\|_2$$
4. **Custo e Massa Final (equação de Tsiolkovsky):**
   $$\Delta V_{\text{total}} = \sum_{k=0}^{M-1} \Delta v_k, \qquad m_f = 1500.0 \exp\!\left(-\frac{\Delta V_{\text{total}}}{g_0 I_{\text{sp}}}\right)$$
5. **Violações e Penalização Proporcional (gradiente contínuo de retorno à região viável):**
   - $v_{\text{launch}} = \max(0,\; v_{\infty,0} - 4000.0)\;\text{(m/s)}$
   - $v_{\text{dur}} = \max(0,\; (t_{M-1} - t_0) - 3652.5)\;\text{(dias)}$
   - $v_{\text{mass}} = \max(0,\; 500.0 - m_f)\;\text{(kg)}$
   - $v_{\text{tof}} = \sum_{a=0}^{M-2} \left[\max(0,\; \text{ToF}_{\min} - \Delta t_a) + \max(0,\; \Delta t_a - \text{ToF}_{\max})\right]\;\text{(dias)}$, onde $\Delta t_a = t_{a+1} - t_a$
   - $n_{\text{fail}} =$ número de arcos com falha numérica em Lambert (`!ok` ou valores não-finitos)
   $$F(X) = \Delta V_{\text{total}} + W_{\text{launch}}\, v_{\text{launch}} + W_{\text{dur}}\, v_{\text{dur}} + W_{\text{mass}}\, v_{\text{mass}} + W_{\text{tof}}\, v_{\text{tof}} + W_{\text{fail}}\, n_{\text{fail}}$$
   Com pesos iniciais recomendados: $W_{\text{launch}} = 10.0$, $W_{\text{dur}} = 500.0$, $W_{\text{mass}} = 200.0$, $W_{\text{tof}} = 1000.0$, $W_{\text{fail}} = 10^6$. Uma solução é marcada `feasible = true` **se e somente se** $n_{\text{fail}} = 0$ e todas as violações ($v_{\text{launch}}, v_{\text{dur}}, v_{\text{mass}}, v_{\text{tof}}$ e limites $[L_k, U_k]$) forem exatamente $0$.

#### Caso 2: GTOC11 ($M \in [36, 43]$)
1. **Nó de Partida ($k = 0$, Terra):**
   $$v_{\infty,0} = \|\mathbf{v}_{\text{dep}}[0] - \mathbf{V}_{\text{Terra}}(t_0)\|_2, \qquad \Delta v_0 = \max(0,\; v_{\infty,0} - 6000.0)$$
   *(Reproduz exatamente `compute_rv_dv` de `DP_flyby_2impulse_user_function.cpp`, onde o excesso além de $6000\;\text{m/s}$ é somado ao $\Delta V$ — ex.: `10.84 m/s` na Missão 1 de `DP_result.txt`; suportar opcionalmente flag `--strict-launch-vinf` para penalizar $v_{\infty,0} > 6000.0$ se desejado em análise complementar).*
2. **Nós de Sobrevoo Intermediários ($k = 1 \dots M-2$):**
   - Calcular $\mathbf{d}_{\text{in}} = \mathbf{V}_{b_k}(t_k) - \mathbf{v}_{\text{arr}}[k-1]$ e $\mathbf{d}_{\text{out}} = \mathbf{v}_{\text{dep}}[k] - \mathbf{V}_{b_k}(t_k)$.
   - Aplicar o método cônico analítico (`conic_flyby_dv`, refatorado de `rendezvous2flyby` com raio $r = (2.0 - 10^{-12})\cdot 1000\;\text{m/s}$), resolvendo $\Delta v_k = \min_{\|\mathbf{v}_{\text{flyby}} - \mathbf{V}_{b_k}\| \le 2000} (\|\mathbf{v}_{\text{arr}}[k-1] - \mathbf{v}_{\text{flyby}}\| + \|\mathbf{v}_{\text{flyby}} - \mathbf{v}_{\text{dep}}[k]\|)$.
3. **Último Alvo ($k = M-1$):**
   $$\Delta v_{M-1} = \max\!\left(0,\; \|\mathbf{V}_{b_{M-1}}(t_{M-1}) - \mathbf{v}_{\text{arr}}[M-2]\|_2 - 2000.0\right)$$
4. **Penalizações:** $v_{\text{tof}}$ ($\text{ToF}_a \in [1.0,\; 730.5]$ dias), $v_{\text{win}} = \sum_{k=0}^{M-1} [\max(0, L_k - t_k) + \max(0, t_k - U_k)]$ (além do reparo que mantém $t_k \in [L_k, U_k]$) e $W_{\text{fail}}\, n_{\text{fail}}$.

---

### 2.2 Cromossomo, Reparo Cronológico e Cache Incremental

1. **Cromossomo:** `std::vector<double> t` de dimensão $M$ (épocas em MJD).
2. **Operador Determinístico de Reparo Cronológico (`repair_chronological`):**
   Garante $L_k \le t_k \le U_k$ e $t_k + \text{ToF}_{\min} \le t_{k+1}$ em tempo $O(M)$ sem perder viabilidade de fronteira:
   - Pré-condição validada na instância: para todo $k$, $U'_{k} = \min(U_k, U'_{k+1} - \text{ToF}_{\min})$ e $L'_{k} = \max(L_k, L'_{k-1} + \text{ToF}_{\min})$ satisfazem $L'_k \le U'_k$.
   - Passagem *forward*: para $k = 0 \dots M-1$, definir `lo = (k == 0) ? L'[0] : std::max(L'[k], t[k - 1] + ToF_min)` e fazer `t[k] = std::clamp(t[k], lo, U'[k])`.
3. **Estrutura de Cache Incremental (`Individual`):**
   - Por corpo ($k = 0 \dots M-1$): `R[k]` (`Vec3`), `V[k]` (`Vec3`), `dv_node[k]` (`double`).
   - Por arco ($a = 0 \dots M-2$): `v_dep[a]` (`Vec3`), `v_arr[a]` (`Vec3`), `arc_ok[a]` (`bool`), `arc_tof[a]` (`double`).
   - Métricas agregadas: `raw_dv`, `final_mass`, `penalty`, `fitness`, `feasible`, `cache_valid`.
4. **Atualização Incremental em Perturbação de 1 Gene ($t_k \to t'_k$):**
   - **1 propagação kepleriana:** recalcular apenas `(R[k], V[k])` no instante $t'_k$.
   - **No máximo 2 chamadas Lambert:** arco $k-1$ (se $k > 0$) e arco $k$ (se $k < M-1$).
   - **Recálculo local de `dv_node` nos até 3 nós acoplados:**
     - Nó $k-1$ (se $k > 0$), pois $\mathbf{v}_{\text{dep}}[k-1]$ mudou;
     - Nó $k$ (sempre), pois $(R[k], V[k])$, $\mathbf{v}_{\text{arr}}[k-1]$ e $\mathbf{v}_{\text{dep}}[k]$ mudaram;
     - Nó $k+1$ (se $k < M-1$), pois $\mathbf{v}_{\text{arr}}[k]$ mudou!
   - **Atualização da soma e penalidades:** recomputar `raw_dv = sum(dv_node)` e `penalty` em $O(M)$ (custo desprezível comparado a Lambert e **elimina erro acumulativo de ponto flutuante**).
   - **Rollback em $O(1)$:** antes de testar $t'_k$ na Busca Local, salvar apenas `t[k]`, `(R[k], V[k])`, até 2 `ArcCache` (`k-1`, `k`), até 3 `dv_node` (`k-1`, `k`, `k+1`) e os escalares agregados. Se o movimento não melhorar, restaurar o snapshot sem nenhuma chamada adicional de Lambert.

---

### 2.3 Operadores do AG Real e Extensão Memética

| Componente | Especificação Padrão | Observações de Implementação |
|---|---|---|
| **População & Inicialização** | $N_{\text{pop}} = 60\text{–}80$ (configurável via CLI) | Mistura de 50% amostragem uniforme em $[L_k, U_k]$ + 50% perfis quase-lineares com ruído gaussiano, todos reparados via `repair_chronological`. |
| **Seleção** | Torneio com $K = 3$ | Compara `fitness`; em caso de empate exato, desempata por `raw_dv` e viabilidade. |
| **Elitismo** | $N_{\text{elite}} = 2$ | Copia indivíduos elite intactos (incluindo cache já avaliado: **0 chamadas Lambert extras**). |
| **Crossover BLX-$\alpha$** | $P_{\text{cross}} = 0{,}85$, $\alpha = 0{,}3$ | Para cada gene $k$: $I_k = |x_k^{(1)} - x_k^{(2)}|$, amostra uniforme em $[\min - \alpha I_k,\, \max + \alpha I_k]$, seguida de `repair_chronological`. |
| **Mutação Gaussiana** | $P_{\text{mut}} = 2{,}0 / M$ por gene | $\sigma(p) = \sigma_{\text{init}} \cdot (\sigma_{\text{final}} / \sigma_{\text{init}})^p$, com $p = \text{LambertsUsados}/\text{Max\_Lamberts} \in [0, 1]$, $\sigma_{\text{init}} = 40{,}0\;\text{d} \to \sigma_{\text{final}} = 0{,}5\;\text{d}$, seguida de `repair_chronological`. |
| **Controle de Estagnação** | Gatilho após $G_{\text{stag}} = 25\text{–}30$ gerações sem melhora | **Imigração parcial:** substitui os 30% piores indivíduos por novos indivíduos gerados e reparados aleatoriamente (respeitando `Max_Lamberts`). |
| **Busca Local Incremental (Núcleo Memético)** | Lamarckiana seletiva a cada 10 gerações + polimento final | **Periódica:** aplicada ao melhor indivíduo + 1 indivíduo aleatório do top 25% (excluindo o melhor), via *Pattern Search* por coordenadas em ordem embaralhada (`std::shuffle`), passos `[8.0, 2.0, 0.5]` dias dentro de $[\max(L_k, t_{k-1}+\text{ToF}_{\min}),\, \min(U_k, t_{k+1}-\text{ToF}_{\min})]$ (com teto de orçamento por chamada, ex. 1 varredura completa por passo).<br>**Polimento final:** aplicado no melhor global ao atingir a reserva final (ex. últimos 10% do orçamento ou passo a passo) com passos `[2.0, 0.5, 0.1, 0.02]` dias. |
| **Orçamento Justo (`Max_Lamberts`)** | Contador global único por execução | 1 unidade = 1 chamada a `solve_lambert_arc`. Avaliação completa = $M-1$ Lamberts; passo incremental = 1 ou 2 Lamberts. AG Puro e AG Memético encerram imediatamente ao atingir `Max_Lamberts`. |

---

### 2.4 Árvore de Diretórios Proposta em `Codigo/`

```text
Codigo/
├── CMakeLists.txt
├── include/
│   ├── astro/
│   │   ├── constants.hpp          # Constantes físicas e fatores de conversão exatos
│   │   ├── orbit_math.hpp         # Vec3, álgebra vetorial sem alocações em heap
│   │   ├── kepler.hpp             # E2f, E2M, f2E, M2E, f0dt2ft, coe2rv, rv2coe
│   │   ├── lambert.hpp            # Izzo/LambEval limpo + wrapper com LambertBudget
│   │   └── conic_flyby.hpp        # Método cônico (rendezvous2flyby sem new/delete)
│   ├── problem/
│   │   ├── catalog.hpp            # Estruturas e leitores GTOC4_data / Candidate_Asteroids
│   │   ├── mission_instance.hpp   # Sequências fixas, épocas da Terra, janelas [L_k, U_k], BKS
│   │   └── evaluator.hpp          # Avaliação completa e atualização incremental + snapshot
│   ├── ga/
│   │   ├── config.hpp             # Parâmetros do AG, Memético e orçamento Max_Lamberts
│   │   ├── individual.hpp         # Cromossomo + Cache por nó e por arco
│   │   ├── repair.hpp             # Operador determinístico de reparo cronológico
│   │   ├── operators.hpp          # Inicialização, Torneio, BLX-alpha, Mutação Gaussiana
│   │   ├── local_search.hpp       # Pattern Search incremental embaralhado
│   │   └── genetic_algorithm.hpp  # Loop evolutivo (AG Puro e AG Memético), elitismo, imigração
│   └── experiment/
│       ├── runner.hpp             # Orquestração de 30 rodadas determinísticas por instância
│       ├── stats.hpp              # Melhor, Média, Mediana, Desvio Padrão, Gap % vs. BKS
│       └── exporter.hpp           # Exportação CSV/TXT (trajetória, estatísticas, convergência)
├── src/
│   ├── astro/                     # Implementações (.cpp) dos módulos de astrodinâmica
│   ├── problem/                   # Parsers, instâncias e avaliador
│   ├── ga/                        # Implementação dos operadores, busca local e AG
│   ├── experiment/                # Runner, estatísticas e exportadores
│   └── main.cpp                   # Interface CLI do executável principal
├── tests/
│   ├── test_sanity_physics.cpp    # GATE: reprodução GTOC4 + 10 missões GTOC11
│   ├── test_incremental_cache.cpp # Equivalência completa vs. incremental (tol <= 1e-9)
│   ├── test_repair_operators.cpp  # Invariantes de reparo, BLX-alpha e mutação
│   └── test_budget_determinism.cpp# Reprodutibilidade por seed e respeito estrito a Max_Lamberts
├── scripts/
│   └── analyze_results.py         # (Opcional) Gráficos de convergência e tabelas comparativas
└── results/                       # Saídas experimentais (gerado em execução)
```

---

## 3. Fases de Desenvolvimento (Milestones)

### Fase 0 — Preparação do Ambiente e Build System (`Codigo/`)
- **Dependências:** nenhuma (início do projeto).
- **Tarefas:**
  - [ ] Renomear/remover a pasta vazia `Código/` (com acento) e criar a estrutura de diretórios `Codigo/`.
  - [ ] Criar `Codigo/CMakeLists.txt` configurando C++20 (compatível com C++17), flags `-O3 -Wall -Wextra` em modo `Release`, biblioteca central `traj_opt_core`, executável principal `gtoc_opt` e integração com `enable_testing()` / CTest.
  - [ ] Validar configuração e build fora da árvore em `Codigo/build/`.
- **Critério de Aceite:** `cmake -S Codigo -B Codigo/build -DCMAKE_BUILD_TYPE=Release && cmake --build Codigo/build -j` conclui com código 0 e **zero warnings**.

### Fase 1 — Núcleo de Astrodinâmica, Parsers e Definição das Instâncias
- **Dependências:** Fase 0 concluída.
- **Tarefas:**
  - [ ] Portar e limpar `constants.hpp`, `orbit_math.hpp`, `kepler.hpp/.cpp`, `lambert.hpp/.cpp` e `conic_flyby.hpp/.cpp` a partir de `Multi_Flyby_Dynamic_Programming/code/include` e `GTOC11_problem.cpp`.
  - [ ] Corrigir no porte: (i) eliminar alocações `new`/`delete` em `rendezvous2flyby`/`calc_transF2X`; (ii) adicionar guarda `iter < MAX_ITER` no `while` de `find_minimum`; (iii) ignorar ramos com `flag0 != 1` ou valores não-finitos em `LambEval`; (iv) corrigir `V_Norm2(r1,3)==0 && V_Norm2(r2,3)==0`.
  - [ ] Implementar `LambertBudget` (contador centralizado dentro da função única de resolução de arco).
  - [ ] Implementar parsers para `GTOC4_data.txt` (pulando 2 linhas `#`, propagando asteroides com época `!= 54800.0` MJD para `54800.0` MJD) e `Candidate_Asteroids.txt` (83.453 linhas TAB).
  - [ ] Codificar sequências fixas e janelas `[L_k, U_k]` para GTOC4 Impulsivo ($\pm 480$ d) e GTOC11 Missões 1–10 ($\pm 600$ d limitado por $T_{\max,k}$).
- **Critério de Aceite:** executável/teste auxiliar carrega exatamente 1.436 asteroides (GTOC4) e 83.453 asteroides (GTOC11), confirma que o asteroide fora da época `54800.0` é propagado corretamente, e verifica que todas as datas de referência (`T_GTOC4`, BKS Tabela 2 e `T_GTOC11_multi` / `DP_result.txt`) estão dentro de `[L_k, U_k]`.

### Fase 2 — [GATE OBRIGATÓRIO] Avaliador Físico Completo e Teste de Sanidade Física
- **Dependências:** Fase 1 concluída.
- **Tarefas:**
  - [ ] Implementar `Evaluator::evaluate_full` para **GTOC4 Impulsivo** (incluindo $v_{\infty,0}$, $\Delta v_k$, $m_f$, duração, limites de ToF e penalidades proporcionais) e para **GTOC11** (incluindo método cônico nos intermediários e $\Delta v_{M-1}$ no alvo final).
  - [ ] Implementar `tests/test_sanity_physics.cpp` e registrá-lo no CTest.
  - [ ] Validar 3 conjuntos de referência:
    1. GTOC4 com `T_GTOC4` (`chimpionship_result`): reproduzir `26.167,447 m/s` e `616,325 kg`.
    2. GTOC4 com as 50 épocas da Tabela 2 de `2508.02904v1.pdf`: confirmar $\Delta V \approx 25.210{,}4\;\text{m/s}$ e $m_f \approx 636{,}7\;\text{kg}$ (diferença residual exclusivamente pelo arredondamento a $0{,}1\;\text{dia}$ na tabela do artigo; validar também convergência à BKS exata após ajuste fino local $\le 0{,}05\;\text{d}$).
    3. GTOC11 Missões 1–10 com as épocas de `DP_result.txt` e com `T_GTOC11_multi`: reproduzir `13.944,823 m/s` (Missão 1) e todas as demais 9 missões com erro `< 0,5 m/s` (e `< 1e-2 m/s` quando usadas épocas de precisão dupla).
- **Critério de Aceite (Gate Bloqueante):** `ctest -R test_sanity_physics --output-on-failure` passa 100%. **É proibido avançar para a Fase 3 antes deste gate passar.**

### Fase 3 — Cromossomo, Reparo Cronológico e Cache Incremental
- **Dependências:** Gate da Fase 2 aprovado.
- **Tarefas:**
  - [ ] Implementar `Individual` com cache por corpo (`R`, `V`, `dv_node`) e por arco (`v_dep`, `v_arr`, `arc_ok`, `arc_tof`).
  - [ ] Implementar `repair_chronological` determinístico em duas passagens ($O(M)$).
  - [ ] Implementar `Evaluator::evaluate_single_date_update` com atualização de 1 efeméride, no máximo 2 Lamberts e até 3 nós (`k-1, k, k+1`), mais mecanismo de snapshot/rollback em $O(1)$.
  - [ ] Criar `tests/test_repair_operators.cpp` e `tests/test_incremental_cache.cpp`.
- **Critério de Aceite:**
  - `test_repair_operators`: 1.000 vetores aleatórios (inclusive invertidos e fora de `[L_k, U_k]`) após reparo satisfazem $L_k \le t_k \le U_k$ e $t_{k+1} - t_k \ge \text{ToF}_{\min}$.
  - `test_incremental_cache`: 1.000 perturbações aleatórias em GTOC4 e GTOC11 (incluindo bordas $k=0$ e $k=M-1$) satisfazem $|\Delta V_{\text{inc}} - \Delta V_{\text{full}}| \le 10^{-9}\;\text{m/s}$ e $|F_{\text{inc}} - F_{\text{full}}| \le 10^{-9}$, consumindo exatamente 1 Lambert nas bordas e 2 no interior, além de rollback restaurar o estado bit-a-bit.

### Fase 4 — Algoritmo Genético Puro (Codificação Real) e Orçamento `Max_Lamberts`
- **Dependências:** Fase 3 concluída.
- **Tarefas:**
  - [ ] Implementar inicialização híbrida, Seleção por Torneio ($K = 3$), Elitismo ($N_{\text{elite}} = 2$), Crossover BLX-$\alpha$ ($P_{\text{cross}} = 0{,}85$, $\alpha = 0{,}3$) e Mutação Gaussiana adaptativa ($P_{\text{mut}} = 2/M$, $\sigma$ decrescente com o orçamento).
  - [ ] Implementar loop evolutivo do **AG Puro** com controle de estagnação e imigração parcial (30% dos piores após 25–30 gerações sem melhora).
  - [ ] Integrar parada rigorosa pelo contador global `Max_Lamberts` e registro de histórico (`LambertsUsados`, `BestFitness`, `BestRawDV`, `Feasible`).
  - [ ] Criar `tests/test_budget_determinism.cpp`.
- **Critério de Aceite:** `test_budget_determinism` verifica: (1) `lamberts_used <= Max_Lamberts`; (2) duas execuções com a mesma seed `std::mt19937_64` produzem trajetórias idênticas bit-a-bit; (3) seeds distintas divergem; (4) `BestFitness` é monotonicamente não-crescente ao longo das gerações.

### Fase 5 — Extensão Memética (Busca Local Incremental)
- **Dependências:** Fase 4 concluída.
- **Tarefas:**
  - [ ] Implementar `local_search_pattern` usando `evaluate_single_date_update` + rollback $O(1)$, ordem de coordenadas embaralhada por `std::shuffle`, respeitando $[\max(L_k, t_{k-1}+\text{ToF}_{\min}),\, \min(U_k, t_{k+1}-\text{ToF}_{\min})]$.
  - [ ] Integrar no loop do **AG Memético**: a cada 10 gerações aplicar no melhor indivíduo + 1 aleatório do top 25% com passos `[8.0, 2.0, 0.5]` dias, e aplicar polimento final no melhor global com passos `[2.0, 0.5, 0.1, 0.02]` dias dentro do mesmo `Max_Lamberts`.
  - [ ] Estender `test_budget_determinism` para cobrir o modo `memetic`.
- **Critério de Aceite:** no modo `memetic`, `lamberts_used <= Max_Lamberts` é respeitado estritamente, a reprodutibilidade por seed é preservada e uma execução piloto em GTOC4 Impulsivo e GTOC11 Missão 1 demonstra melhora de convergência/qualidade frente ao AG Puro sob o mesmo orçamento.

### Fase 6 — CLI, Runner Experimental (30 Seeds), Exportação e Estudo de Ablação
- **Dependências:** Fase 5 concluída.
- **Tarefas:**
  - [ ] Implementar CLI completa em `src/main.cpp`: `--problem gtoc4|gtoc11`, `--mission 1..10|all`, `--mode ga|memetic|both`, `--budget <N>`, `--runs 30`, `--seed-base <S>`, `--pop <N>`, `--out-dir <path>`.
  - [ ] Implementar `Runner` para executar 30 rodadas independentes com sementes determinísticas `seed = seed_base + run_id` (suportando paralelismo opcional entre rodadas independentes via OpenMP mantendo o RNG isolado por rodada).
  - [ ] Implementar `Exporter` gerando por instância e modo:
    - `summary_stats.csv` (Melhor, Média, Mediana, Desvio Padrão, Gap % vs. BKS, Taxa de Viabilidade, Tempo Médio s);
    - `runs_detail.csv` (métrica final, `LambertsUsados`, tempo e seed de cada uma das 30 rodadas);
    - `best_trajectory.csv` e `best_trajectory.txt` (formato legível compatível com `DP_result.txt`);
    - `convergence_run_<id>.csv` (`LambertsUsados`, `BestFitness`, `BestRawDV`, `Feasible`).
  - [ ] Executar campanha experimental completa (AG Puro vs. AG Memético sob mesmo `Max_Lamberts`) e gerar resumo comparativo em `Codigo/results/`.
- **Critério de Aceite:** toda a suíte `ctest` passa; execução de 30 seeds gera todos os arquivos CSV/TXT íntegros sem `NaN`, com comparação justa sob idêntico `Max_Lamberts` entre AG Puro e AG Memético.

---

## 4. Riscos Técnicos e Mitigação

| Risco | Causa Raiz Observada no Código de Referência | Mitigação Obrigatória em `Codigo/` |
|---|---|---|
| **Falhas e singularidades do solucionador de Lambert** | Em `Lambert.cpp`: `LambEval` tem `if (flag0 != 1 || e >= 1.0) continue;` comentado (linha 344) e retorna velocidades artificiais `10000.0` se `iter >= Maxiter`; transferências com $\theta \approx 0$ ou $\pi$ ou $t_f \le 0$ podem gerar `NaN` ou velocidades espúrias. | Exigir `flag0 == 1` e `std::isfinite` em `v0, v1`; manter $\Delta t_a \ge \text{ToF}_{\min}$ via reparo; atribuir `arc_ok = false` e penalidade $W_{\text{fail}}$ quando nenhum ramo convergir; calcular $N_{\max}$ com guarda para evitar overflow inteiro. |
| **Divergência entre avaliação completa e incremental** | Perturbar $t_k$ modifica o arco $k-1$ e o arco $k$, afetando o $\Delta v$ de **três nós** (`k-1`, `k`, `k+1`) — especialmente crítico no método cônico do GTOC11. Somas incrementais `total += new - old` acumulam drift de ponto flutuante. | Atualizar os até 3 nós afetados (`k-1, k, k+1`), recomputar `raw_dv = sum(dv_node)` em $O(M)$, usar snapshot completo das partes alteradas no rollback e validar na suíte `test_incremental_cache` com tolerância $\le 10^{-9}$. |
| **Bugs herdados nos arquivos de referência** | `find_minimum` (em `GTOC11_problem.cpp`) usa `while (fabs(fc - fd) > TOLERANCE)` sem verificar `iter < MAX_ITER`; `rendezvous2flyby`/`calc_transF2X` fazem `new`/`delete` em toda chamada; linha 389 testa `V_Norm2(r1,3)==0` duas vezes; Terra do GTOC4 já traz anomalia verdadeira `f` e não `M`; 1 asteroide do GTOC4 está na época `54745.918153`. | Corrigir todos esses pontos no porte para `Codigo/include/astro` e `Codigo/src/problem`, substituindo matrizes alocadas por `std::array<std::array<double, 3>, 3>` e verificando o gate da Fase 2. |
| **Vazamento de orçamento ou perda de determinismo** | Uso de `rand()`, RNG global compartilhado ou chamadas a Lambert fora do contador (ex.: durante imigração ou busca local). | Passar `std::mt19937_64&` explicitamente; contabilizar `budget.consume()` dentro da única função `solve_lambert_arc()`; testar bit-a-bit em `test_budget_determinism`. |
| **Custo computacional nas 30 rodadas × 11 instâncias** | Avaliações excessivas com muitas multi-revoluções e alocações repetidas de memória em `std::vector`. | Pré-alocar buffers do `Individual` uma única vez, usar `std::array<double, 3>`, compilar com `-O3`, usar cache incremental na busca local (máx. 2 Lamberts/teste) e permitir paralelismo exclusivamente *entre* rodadas independentes. |
