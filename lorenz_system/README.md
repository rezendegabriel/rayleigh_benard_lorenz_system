# lorenz_system



# Simulações do Sistema de Lorenz em Julia

Este repositório contém três códigos em Julia utilizados para a realização 
das simulações numéricas do Sistema de Lorenz, empregando o método de 
Runge–Kutta de quarta ordem (RK4).

Os scripts geram automaticamente todas as figuras utilizadas no trabalho, 
incluindo trajetórias tridimensionais do atrator, projeções bidimensionais 
do espaço de fases, análise da sensibilidade às condições iniciais e 
comparações temporais das variáveis do modelo.

## Arquivos

- **SdL.jl**  
  Realiza a simulação do sistema de Lorenz para **seis condições iniciais distintas**, 
  gerando os atratores em 3D e as projeções nos planos \(x \times y\), 
  \(x \times z\) e \(y \times z\).

- **SdL_ST.jl**  
  Gera as **séries temporais das variáveis** \(x(t)\), \(y(t)\) e \(z(t)\),
  comparando a trajetória sem perturbação com uma trajetória ligeiramente perturbada.

- **SdLTT.jl**  
  Realiza o estudo da **sensibilidade às condições iniciais**, introduzindo 
  pequenas perturbações nas coordenadas \(x\), \(y\) e \(z\), e gerando 
  figuras tridimensionais comparativas entre a trajetória base e as trajetórias perturbadas.

## Execução

Para executar os scripts, é necessário ter a linguagem Julia instalada e o pacote `Plots` disponível.


