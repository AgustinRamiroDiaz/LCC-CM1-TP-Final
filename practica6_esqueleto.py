#! /usr/bin/python

# 6ta Practica Laboratorio 
# Complementos Matematicos I
# Ejemplo parseo argumentos

import argparse

import matplotlib.pyplot as plt
import numpy as np


class LayoutGraph:

    def __init__(self, grafo, iters, refresh, c1, c2, verbose = False):
        """Parametros de layout:
        iters: cantidad de iteraciones a realizar
        refresh: Numero de iteraciones entre actualizaciones de pantalla.
        0 -> se grafica solo al final.
        c1: constante usada para calcular la repulsion entre nodos
        c2: constante usada para calcular la atraccion de aristas"""

        # Guardo el grafo
        self.grafo = grafo
        self.nodos = self.grafo[0]
        self.aristas = self.grafo[1]

        # Inicializo estado
        # Completar
        self.posiciones = {}
        self.fuerzas = {}
        self.accumx = {}
        self.accumy = {}

        # Guardo opciones
        self.iters = int(iters)
        self.verbose = verbose
        # TODO: faltan opciones
        self.refresh = refresh
        self.c1 = c1
        self.c2 = c2
        self.c = 1
        self.area = 100
        self.k = self.c * np.sqrt(self.area / len(self.nodos))

        pass

    def layout(self):
        """
        Aplica el algoritmo de Fruchtermann-Reingold para obtener (y mostrar)
        un layout
        """
        self.algoritmoFruchtermanReingold()
        pass

    def posicionesAleatorias(self):
        for vertice in self.grafo[0]:
            self.posiciones[vertice] = [np.random.random_sample(), np.random.random_sample()]
        pass

    def plotear(self):
        plt.clf()
        for ni, nj in self.aristas:
            posicionOrigen = self.posiciones[ni]
            posicionDestino = self.posiciones[nj]

            absisaOrigen = posicionOrigen[0]
            ordenadaOrigen = posicionOrigen[1]
            absisaDestino = posicionDestino[0]
            ordenadaDestino = posicionDestino[1]
            plt.plot([absisaOrigen, absisaDestino], [ordenadaOrigen, ordenadaDestino])

        plt.pause(0.01)
        pass

    def algoritmoFruchtermanReingold(self):
        # Seteamos posiciones iniciales aleatorias
        self.posicionesAleatorias()

        accum = {}
        for k in range(self.iters):
            self.step()
        plt.show()
        pass

    def step(self):
        self.initializeAccumulators()
        self.computeAttractionForces()
        self.computeRepulsionForces()
        self.updatePositions()
        self.plotear()
        pass

    def initializeAccumulators(self):
        for vertice in self.grafo[0]:
            self.accumx[vertice] = 0
        for vertice in self.grafo[0]:
            self.accumy[vertice] = 0
        pass

    def computeAttractionForces(self):
        for ni, nj in self.aristas:
            distance = distanciaEuclidiana(self.posiciones[ni], self.posiciones[nj])
            modfa = self.attraction(distance)
            fx = modfa * (self.abcisa(nj) - self.abcisa(ni)) / distance
            fy = modfa * (self.ordenada(nj) - self.ordenada(ni)) / distance

            self.accumx[ni] += fx
            self.accumy[ni] += fy

            self.accumx[nj] -= fx
            self.accumy[nj] -= fy
        pass

    def computeRepulsionForces(self):
        for ni in self.nodos:
            for nj in self.nodos:
                if ni != nj:
                    distance = distanciaEuclidiana(self.posiciones[ni], self.posiciones[nj])
                    modfa = self.attraction(distance)
                    fx = modfa * (self.abcisa(nj) - self.abcisa(ni)) / distance
                    fy = modfa * (self.ordenada(nj) - self.ordenada(ni)) / distance

                    self.accumx[ni] += fx
                    self.accumy[ni] += fy

                    self.accumx[nj] -= fx
                    self.accumy[nj] -= fy
        pass

    def updatePositions(self):
        # TODO: revisar bordes ventana (?)
        for node in self.nodos:
            self.posiciones[node][0] += self.accumx[node]
            self.posiciones[node][1] += self.accumy[node]
        pass

    def attraction(self, distance):
        return distance ** 2 / self.k

    def repultion(self, distance):
        return self.k ** 2 / distance

    def abcisa(self, v):
        return self.posiciones[v][0]

    def ordenada(self, v):
        return self.posiciones[v][1]


def distanciaEuclidiana(ni, nj):
    return np.sqrt((ni[0] - nj[0]) ** 2 + (ni[1] - nj[1]) ** 2)


assert distanciaEuclidiana((0, 0), (3, 4)) == 5


def leeGrafoArchivo(file_path):
    archivo = open(file_path)
    lista = archivo.readlines()
    cantidadDeVertices = int(lista.pop(0))
    listaVertices = []
    for x in range(cantidadDeVertices):
        listaVertices.append(lista[x].rstrip("\n"))
    ARITAS = lista[cantidadDeVertices:]
    aristas = []
    for par in ARITAS:
        aristas.append((par.split()[0], par.split()[1]))
    return listaVertices, aristas


def main():
    # Definimos los argumentos de linea de comando que aceptamos
    parser = argparse.ArgumentParser()

    # Verbosidad, opcional, False por defecto
    parser.add_argument(
        '-v', '--verbose',
        action = 'store_true',
        help = 'Muestra mas informacion al correr el programa'
    )
    # Cantidad de iteraciones, opcional, 50 por defecto
    parser.add_argument(
        '--iters',
        type = int,
        help = 'Cantidad de iteraciones a efectuar',
        default = 50
    )
    # Temperatura inicial
    parser.add_argument(
        '--temp',
        type = float,
        help = 'Temperatura inicial',
        default = 100.0
    )
    # Archivo del cual leer el grafo
    parser.add_argument(
        'file_name',
        help = "Archivo del cual leer el grafo a dibujar"
    )

    parser.add_argument(
        'iters',
        help = "Cantidad de iteraciones del algoritmo"
    )

    args = parser.parse_args()

    # Descomentar abajo para ver funcionamiento de argparse
    # print args.verbose
    # print args.iters    
    # print args.file_name
    # print args.temp
    # return

    grafo = leeGrafoArchivo(args.file_name)

    # Creamos nuestro objeto LayoutGraph
    layout_gr = LayoutGraph(
        grafo,  # TODO: Cambiar para usar grafo leido de archivo
        iters = args.iters,
        refresh = 1,
        c1 = 0.1,
        c2 = 5.0,
        verbose = args.verbose
    )

    # Ejecutamos el layout
    layout_gr.layout()
    return


if __name__ == '__main__':
    main()
