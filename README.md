# Introdução

Esse projeto consiste na aplicação da Análise de Componentes Principais à comparação de sequências de nucleotídeos do gene Citocromo C1 (COI), que está presente no DNA mitocondrial de diferentes espécies de animais.

O DNA funciona como o "código de barras" da vida, pois cada parte dele codifica um gene, a unidade fundamental da hereditariedade e que é responsável por determinar a expressão de uma característica biológica. Todos os seres vivos possuem uma molécula de DNA dentro de cada uma de suas células, e essa molécula tem o formato de uma dupla-hélice, cada uma formada por uma sequência de **nucleotídeos**, que por sua vez são formados por uma pentose (uma molécula orgânica), fosfato (um íon) e uma base nitrogenada que pode ser Adenina (A), Timina (T), Citosina (C) ou guanina (G). É comum identificarmos cada nucleotídeo que compõe uma hélice do DNA através de sua base nitrogenada, e é assim que representamos amostras dessa molécula geralmente, como uma _string_ formada por caracteres A, T, G e C. 

A mitocôndria é uma importante organela presente em células eucariontes (que possuem um núcleo envolvendo seu material genético) que, diferente das outras, possui DNA próprio, o que a torna interessante para comparações evolutivas. A escolha do do gene COI para esse trabalho se deve ao fato dele estar presente em um variado grupo de espécies do reino animal, de modo a facilitar a busca de amostras para usarmos aqui.

O intuito desse projeto é comparar sequências do gene COI de diferentes animais e comparar, a luz desse gene, o quão duas espécies podem estar próximas evolutivamente, de forma a comparar também com a nossa própria noção usual de proximidade das espécies: "será que o chimpanzé ficará próximo dos seres humanos?" ou "uma mosca ficará próxima de uma abelha?".

Para isso, iremos aplicar a Análise de Componentes Principais (PCA), que é amplamente usado no campo da bioinformática para encontrar relação entre amostras de dados genéticos. O PCA é capaz de receber dados de grande dimensionalidade (neste caso, o tamanho das cadeias de nucleotídeos que compõe o gene) e reduzir ela para poucas dimensões que concentram as maiores variâncias dos dados. Desse modo, tentaremos aplicar o PCA e, de modo visual, observar como o grau de parentesco entre as espécies se apresenta na proximidade dos pontos que as representam sobre uma base ortonormal dada pelas primeiras componentes principais.

# Obtenção e Tratamento dos Dados

As amostras do gene COI foram obtidas do site GeBank do _National Canter for Biotechnology Information_ (NCBI). Cada amostra é uma sequência de caracteres A, T, C e G e a obtemos para diferentes tipos de espécies de animais. Abaixo, estão as espécies selecionadas e tentamos colocar espécies de grupos taxnonômicos próximos (como Ser Humano, Chimpazé e Camundongo como mamíferos, Mosca, Abelha e Caranguejo como Artrópodes, etc.) a fim de comparar os resultados obtidos pelo PCA e observar se espécies de grupos taxonômicos parecidos estão próximas umas das outras.

| Nome Genérico         | Nome Científico         | Grupo     |
| --------------------- | ----------------------- | --------- |
| Ser humano            | Homo sapiens            | Mamífero  |
| Camundongo            | Mus musculus            | Mamífero  |
| Chimpanzé             | Pan paniscus            | Mamífero  |
| Mosca de Fruta        | Drosophila melanogaster | Artrópode |
| Abelha                | Apis mellifera          | Artrópode |
| Caranguejo            | Limulus polyphemus      | Artrópode |
| Peixe zebra           | Danio rerio             | Peixe     |
| Carpa                 | Cyprinus carpio         | Peixe     |
| Danio pérola          | Danio_albolineatus      | Peixe     |
| Sanguessuga           | Hirudo medicinalis      | Anelídeo  |
| Minhoca               | Pontoscolex corethrurus | Anelídeo  |
| Minhoca da terra      | Lumbricus terrestris    | Anelídeo  |
| Avestruz              | Struthio camelus        | Ave       |
| Galinha               | Gallus gallus           | Ave       |
| Falcão                | Falco peregrinus        | Ave       |
| Rã                    | Xenopus laevis          | Anfíbio   |
| Axalote               | Ambystoma mexicanum     | Anfíbio   |
| Salamandra            | Pleurodeles waltl       | Anfíbio   |
| Lesma do mar          | Aplysia californica     | Molusco   |
| Polvo                 | Octopus vulgaris        | Molusco   |
| Caramujo de água doce | Lymnaea stagnalis       | Molusco   |

Após obter as sequências dos genes para as espécies listadas, precisamos tratá-las para um formato numérico que onde o PCA pode ser aplicado. A solução para isso é simples, apenas iteramos sobre cada sequência e associamos um valor para cada base, de modo que nosso dicionário será:

| Base | Valor |
| ---- | ----- |
| A    | 1     |
| C    | 2     |
| G    | 3     |
| T    | 4     |

E agora temos uma estrutura de dados na forma da uma matriz em que cada linha é uma espécie e cada coluna é a posição de um nucleotídeo do gene COI. Com isso, estamos aptos a aplicar o PCA sobre os dados.

# Obtendo as Componentes Principais

Suponha que a matriz obtida seja $M$ tal que $X = \{x_1, x_2, \ldots, x_n\}$ seja o conjunto das espécies. Cada $x_i \in X$ é um vetor que compõe as linhas de $M$. Nosso objetivo é encontrar uma base ortonormal para um subespaço sobre o qual iremos projetar os vetores de $X$ de modo que as normas das projeções sejam tão próximas das normas originais o quanto possível. O PCA será usado justamente para alcançar esse objetivo, de modo que os vetores da base ortonormal que estamos buscando serão os autovetores associados aos maiores autovalores da matriz de covariância $M^t M$, que é uma matriz simétrica, o que nos garante a existência de uma base ortonormal de autovetores. 

A forma escolhida para obter as componentes principais foi através da Decomposição de Valores Singulares (SVD) da nossa matriz $M$ de espécies por nucleotídeos. Realizando a fatoração SVD da matriz $M$, nós obteremos $$ M = E \times \Sigma \times N^t $$em que $E$ e $N$ são matrizes ortogonais e $\Sigma$ é uma matriz diagonal com os valores singulares de $M$. A nossa base ortonormal que gera o subespaço que maximiza as normas dos vetores de $X$ são as linhas de $N^t$, de modo que as $k$ primeiras componentes principais das espécies são as $k$ primeiras linhas de $N^t$. Se $L = \{n_1, \ldots, n_k\}$ é o conjunto das $k$ primeiras principais, tal que para cada $n_j \in L$ temos que $n_j$ é um vetor que compõe umas das $k$ primeiras linhas de $N^t$, então cada entrada $e_{i, j}$ da matriz $E \times \Sigma$ é o tamanho da projeção de $x_i \in X$ em $n_j$, tal que $$ e_{i,j} = x_i^tn_j $$
Portanto, para qualquer dimensão $k$, se quisermos saber o tamanho das projeções dos vetores de $X$ de espécies, basta tomarmos as primeiras $k$ colunas da matriz $E \times \Sigma$. 
# Resultados
Abaixo, está nossa matriz $M$ reconstruída com diferentes quantidades de componentes principais. Aqui se justifica nossa escolha de organizarmos as amostras de sequências de modo que grupos taxonômicos semelhantes estejam em linhas consecutivas. Percebamos, na reconstrução com 2 componentes principais, por exemplo, como que espécies taxonomicamente próximas, salvo algumas exceções, são semelhantes quanto o aspecto das "barras" na direção das colunas, como se cada grupo tivesse um "código de barras". Conforme, aumentamos o número de componentes principais, mais próximos estamos do aspecto da matriz original e mais diferentes vão ficando as espécies individualmente.

![image](https://github.com/user-attachments/assets/abb9d25b-1240-433b-ade8-19d299ecf198)

Para as duas primeiras componentes, percebemos que os mamíferos (que são as primeiras linhas da matriz original) tem bastante similaridade entre si, assim como os artrópodes (que são as linhas seguintes aos mamíferos) e os peixes (que seguem os artrópodes). Já outros grupos, como os moluscos (as últimas linhas da matriz original) ficaram menos parecidos entre si, assim como os anelídeos (que estão abaixo dos peixes).

Na próxima imagem, encontram-se as espécies (com seus nomes genéricos) representadas em duas dimensões cuja base é formada pelas duas primeiras componentes principais.

![image](https://github.com/user-attachments/assets/3779ad25-127d-464f-83e1-a2bdca31d156)


Daqui podemos fazer várias observações, como
- Os mamíferos Ser Humanos e Chimpanzé estão mais próximos entre si do que dos camundongos.

- Os artrópodes mosca, abelha, caranguejo estão bem próximos entre si.

- Os moluscos Polvo e Lesma do mar estão bem róximos quanto a segunda componente, porém mais distantes quanto a primeira. O caramujo de água doce, também um molusco, parece estar distante dos outros nas duas componentes.

- Os anfíbios rã e axalote estão bem próximas tanto na primeira quanto na segunda componente, porém estão da salamandra quanto a segunda componente.

- O peixe-zebra, a carpa e o danio pérola estão mais próximos entre si na segunda componente do que na primeira.

- Os anelídeos minhoca, minhoca da terra e sanguessuga estão mais próximos entre si na primeira componente do que na segunda.

- As aves galinha e avestruz estão mais próximas entre si do que esses do falcão.

Também temos que animais de grupos taxonômicos mais distantes conservam certa semelhança em algumas componentes, como peixe-zebra e avestruz nos PC1 e PC2, sanguessuga e rã no PC2, etc.

Abaixo, representamos as espécies agora em dimensão 3. Podemos observar mais atentamente algumas diferenças a mais nos grupos taxonômicos, como por exemplo os artrópodes abelha e mosca estarem mais perto entre si do que do caranguejo, o que faz sentido, já que os dois primeiros são insetos, enquanto o terceiro é um crustáceo, ou uma galinha e um avestruz serem mais parecidos entre si do que esses de um falcão, o que de fato ocorre se analisarmos um cladograma das aves.

![image](https://github.com/user-attachments/assets/0ce008de-4d2e-43e3-b902-4f78dc1817b3)


No gráfico abaixo, percebemos que boa parte da variância (cerca de 83%) nos dados encontram-se na primeira componente principal. O fit de ajuste do PCA encontra-se nas duas primeiras componentes principais, o que sugere que as variâncias explicadas nas demais componentes são ruídos estatísticos.

![image](https://github.com/user-attachments/assets/0bfadf65-2c61-436e-b349-8dcc97d228de)


| PC  | Variância  |
| --- | ---------- |
| 1   | 0.83124784 |
| 2   | 0.02302604 |
| 3   | 0.01328628 |
| 4   | 0.01207710 |
| 5   | 0.01187563 |
Também é interessante analisar a variância em cada grupo taxonômico, para avaliarmos o quanto as espécies de um mesmo grupo distinguem-se umas das outras. Para isso, iremos observar somente a variância estatística nas duas primeiras componentes principais por esse ser o nosso melhor ajuste.

Abaixo está um gráfico de barras para a variância nos grupos na primeira componente principal.

![image](https://github.com/user-attachments/assets/7a51af2a-c1d4-4da0-83c4-98f0f3bc4040)


Observamos que, por exemplo, os anelídeos e os anfíbios variam muito pouco nessa componente. Isso implica que para as espécies que selecionamos, tanto de anelídeos quanto de anfíbios, elas estão bem próximas quanto a principal característica relacionada ao gene COI. Isso já não vale para os peixes. Para a segunda componente, temos os seguinte gráfico:

![image](https://github.com/user-attachments/assets/9e76f8ad-18f9-41a9-adc8-2baf3dd553ae)


Já para segunda principal característica relacionada ao gene, vemos que os anelídeos variam bem mais entre si dessa vez quando comparados aos demais grupos. Vemos também que para essa segunda característica, mamíferos, artrópodes e peixes, por exemplo, estão bem próximos entre si.

# Conclusão

A análise de componentes principais é uma tecnologia extremamente útil para a bioinformática, sendo amplamente usada para a análise de dados genéticos com os mais diversos fins: análise evolutiva, de expressão gênica, para caracterizar amostras relacionadas a uma determinada condição médica, etc. Aqui vimos de forma simplificada uma aplicação do PCA para comparação de um gene entre várias espécies com o fim de observar como poderíamos interpretar o grau de proximidade evolutiva entre elas com relação a esse gene.

A análise de relações evolutivas é limitada, pois apenas utilizamos um único gene para fazer a análise, além de um número relativamente pequeno de espécies para compará-lo. O gene COI foi escolhido, pois ele está sequenciado em bancos de dados biológicos para várias tipos de espécies, aumentando o sucesso de sua busca para uma espécie qualquer. Inicialmente o projeto tinha a intenção de comparar uma variedade maior de seres, incluindo vegetais, bactérias, vírus, etc. Porém, devido a dificuldade de encontrar um gene comum que esteja sequenciado para vários tipos de organismos, decidiu-se por limitar as comparações somente aos animais. 
