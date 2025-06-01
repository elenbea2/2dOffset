# 2D Offset

Este projeto implementa uma técnica de **offset de polígonos 2D** utilizando C++ e a biblioteca [CGAL](https://www.cgal.org/). 

Este trabalho está sendo desenvolvido durante a minha iniciação científica, no curso de Engenharia de Computação na UTFPR, com orientação dos professores Ricardo Dutra da Silva, Rodrigo Minetto e Neri Volpato.

## 🧱 Tecnologias utilizadas

- C++
- [CGAL](https://www.cgal.org/) (Computational Geometry Algorithms Library)
- CMake
- GCC/G++ (Linux)

# 🛠️ Como compilar

Antes de tudo, certifique-se de que você tem as dependências instaladas:

```bash
sudo apt update
sudo apt install cmake g++ libboost-all-dev libcgal-dev
```

## Clone o repositório
git clone git@github.com:elenbea2/2dOffset.git
cd 2dOffset

## Crie uma pasta de build separada
mkdir build
cd build

## Gere os arquivos de build com CMake
cmake ..

## Compile o projeto
make

## Rode o executável
./offset
