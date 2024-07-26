# Quantum Teleportation Simulation

## Descrição

Este projeto é uma simulação de teletransporte quântico usando o pacote QuTiP (Quantum Toolbox in Python). Ele implementa dois protocolos principais de teletransporte quântico: o Protocolo de Teletransporte Padrão (STP) e o Protocolo de Teletransporte Realista (RTP). O código também inclui a aplicação de diversos modelos de ruído e decoerência para avaliar o desempenho dos protocolos.

## Funcionalidades

- **Protocolo de Teletransporte Padrão (STP):** Simula o teletransporte quântico sem considerar o ruído e a decoerência.
- **Protocolo de Teletransporte Realista (RTP):** Simula o teletransporte quântico levando em conta modelos de ruído e decoerência realistas.
- **Modelos de Ruído:** Inclui canais de Bit-Flip, Phase-Flip, Bit-Phase-Flip, Depolarizing e Amplitude Damping.
- **Modelos de Decoerência:** Inclui o canal de Decoerência de Fase.
- **Análise:** Gera um histograma da fidelidade dos estados quânticos após o teletransporte realista.

## Pré-requisitos

Para executar este projeto, você precisa ter o Python instalado, bem como os seguintes pacotes:

- `qutip`
- `numpy`
- `matplotlib`

Você pode instalar as dependências necessárias com o seguinte comando:

```bash
pip install qutip numpy matplotlib
