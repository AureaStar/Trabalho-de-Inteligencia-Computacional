# Guia Operacional para Agentes de Desenvolvimento (`AGENTS.md`)

Este documento define as regras obrigatórias de trabalho, arquitetura, build, testes e conduta para todos os agentes que implementarem as fases descritas em [ROTEIRO_DESENVOLVIMENTO.md](file:///home/aureacp/Documentos/UFJF/IC/ROTEIRO_DESENVOLVIMENTO.md).

---

## 1. Regras Invioláveis de Workspace

1. **Diretórios estritamente somente leitura:**
   - `Fontes/`
   - `Multi_Flyby_Dynamic_Programming/`
   - **Proibido** criar, editar, renomear, compilar ou gravar logs/saídas dentro dessas duas pastas. Use-as apenas como referência de leitura e fonte dos arquivos de entrada (`GTOC4_data.txt`, `Candidate_Asteroids.txt`).
2. **Diretório único de implementação (`Codigo/`):**
   - Todo código C++, `CMakeLists.txt`, build (`Codigo/build/`), testes, scripts e resultados (`Codigo/results/`) deve ficar exclusivamente dentro de `Codigo/`.
   - Nota sobre a pasta vazia `Código/` (com acento): na Fase 0, padronize para `Codigo/` (sem acento) e nunca misture caminhos acentuados no CMake.

---

## 2. Fluxo de Trabalho por Fases

1. Leia [ROTEIRO_DESENVOLVIMENTO.md](file:///home/aureacp/Documentos/UFJF/IC/ROTEIRO_DESENVOLVIMENTO.md) antes de iniciar qualquer etapa.
2. Execute **uma fase por vez**, seguindo rigorosamente a ordem **Fase 0 $\to$ Fase 6**.
3. Só marque `[x]` no checklist de uma fase em `ROTEIRO_DESENVOLVIMENTO.md` após comprovar todos os critérios de aceite daquela fase por meio de execução real de build e testes.
4. **GATE BLOQUEANTE — Fase 2 (Teste de Sanidade Física):**
   - **É proibido** iniciar a Fase 3 (cache incremental) ou implementar qualquer operador do Algoritmo Genético antes de `test_sanity_physics` passar 100%, reproduzindo as trajetórias e BKS de referência do **GTOC4 Impulsivo** e das **10 missões do GTOC11**.

---

## 3. Padrões de Código e Build

- **Linguagem e Compilador:** C++20 (compatível com C++17), codificação UTF-8, portável para Linux (sem `system("pause")`, sem headers exclusivos do Windows, sem acentos em nomes de arquivos/caminhos).
- **Sistema de Build:** CMake fora da árvore em `Codigo/build/` com flags obrigatórias:
  ```bash
  cmake -S Codigo -B Codigo/build -DCMAKE_BUILD_TYPE=Release
  cmake --build Codigo/build -j
  ```
  O `CMakeLists.txt` deve habilitar `-O3 -Wall -Wextra` no alvo principal e nos testes; o código autoral deve compilar com **zero warnings**.
- **Separação Modular:** mantenha fronteiras estritas entre:
  - `astro/`: física orbital pura (`kepler`, `lambert`, `conic_flyby`, `constants`, `orbit_math`);
  - `problem/`: leitores de catálogo, definição de instâncias/janelas e `Evaluator`;
  - `ga/`: cromossomo/cache, reparo cronológico, operadores, busca local e loop evolutivo;
  - `experiment/`: runner das 30 execuções, estatísticas e exportação CSV/TXT.
- **Boas práticas obrigatórias:**
  - Nunca use `using namespace std;` em arquivos de cabeçalho (`.hpp`).
  - Zero alocações dinâmicas (`new`/`delete`) em rotinas internas de avaliação (ex.: `rendezvous2flyby`, `lambert`, cálculo de $\Delta V$ e busca local).
  - Use exclusivamente `std::mt19937_64` com seeds explícitas passadas por referência. **Proibido** usar `rand()`, `srand()` ou estado RNG global mutável.

---

## 4. Regras de Validação, Testes e Orçamento Justo

1. **Equivalência Numérica do Cache Incremental (Fase 3):**
   - `test_incremental_cache` deve testar perturbações aleatórias em todas as regiões (bordas $k = 0$, $k = M-1$ e interior $0 < k < M-1$) para GTOC4 e GTOC11, exigindo erro máximo:
     $$|\Delta V_{\text{inc}} - \Delta V_{\text{full}}| \le 10^{-9}\;\text{m/s}, \qquad |F_{\text{inc}} - F_{\text{full}}| \le 10^{-9}$$
   - Confirmar contagem de chamadas Lambert: 1 chamada nas bordas e 2 no interior.
2. **Contabilização Exata de `Max_Lamberts` (Fases 4–6):**
   - Toda chamada ao solucionador de arco de Lambert deve passar por um único ponto instrumentado (`LambertBudget`).
   - Tanto **AG Puro** quanto **AG Memético** (incluindo imigração, buscas locais periódicas e polimento final) devem compartilhar a mesma métrica de orçamento e interromper imediatamente ao atingir `Max_Lamberts`.
3. **Reprodutibilidade:**
   - `test_budget_determinism` deve garantir que a mesma seed produza resultado idêntico bit-a-bit nos modos `ga` e `memetic`.
4. **Execução da Suíte de Testes:**
   - Após qualquer alteração em `Codigo/`, rodar:
     ```bash
     ctest --test-dir Codigo/build --output-on-failure
     ```

---

## 5. Restrições do Trabalho Acadêmico (DCC136 — UFJF 2026.3)

- A metaheurística (AG de codificação real, operadores genéticos, reparo cronológico, controle de estagnação e busca local memética) deve ser **integralmente implementada pelo projeto** — é estritamente proibido delegar a otimização para bibliotecas prontas de metaheurísticas ou otimização não-linear externas.
- Apenas as rotinas de física orbital (propagação kepleriana, solucionador de Lambert e método cônico) podem ser adaptadas de `Multi_Flyby_Dynamic_Programming/`, preservando a atribuição de autoria nos cabeçalhos/docstrings.
- O código deve ser limpo, bem estruturado e comentado no *porquê* das decisões, servindo de base direta para a documentação metodológica (lembrando que o relatório final da disciplina não deve conter trechos de código-fonte).

---

## 6. Checklist Rápido ao Concluir Cada Fase

1. Compilar `Codigo/build` em `Release` (`-O3 -Wall -Wextra`) sem warnings.
2. Executar `ctest --test-dir Codigo/build --output-on-failure` com 100% de aprovação.
3. Verificar que nenhum arquivo foi criado ou alterado em `Fontes/` ou `Multi_Flyby_Dynamic_Programming/`.
4. Atualizar os checkboxes da fase concluída em `ROTEIRO_DESENVOLVIMENTO.md` e reportar ao usuário os números e saídas reais observadas nos testes.
