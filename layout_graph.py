#! /usr/bin/python

# 6ta Practica Laboratorio 
# Complementos Matematicos I
# Ejemplo parseo argumentos

import matplotlib.pyplot as plt
import numpy as np
import time
import pprint

class LayoutGraph:

    def __init__(self, grafo, refresh, temperaturaInicial, constanteTemperatura, gravedad, repultionConstant, attractionConstant, iters = 100, verbose = False):
        """Parametros de layout:
        iters: cantidad de iteraciones a realizar
        refresh: Numero de iteraciones entre actualizaciones de pantalla.
        0 -> se grafica solo al final.
        repultionConstant: constante usada para calcular la repulsion entre vértices
        attractionConstant: constante usada para calcular la atraccion de aristas"""

        # Guardo el grafo
        self.grafo = grafo
        self.vertices = self.grafo[0]
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
        if (refresh != 0):
            self.refresh = refresh
        else: 
            self.refresh = self.iters
        self.repultionConstant = repultionConstant
        self.attractionConstant = attractionConstant
        self.c = 3
        self.frameSize = 100
        self.area = self.frameSize ** 2
        self.k = self.c * np.sqrt(self.area / len(self.vertices))
        self.gravity = gravedad
        self.temperaturaInicial = temperaturaInicial
        self.temperatura = self.temperaturaInicial
        self.constanteTemperatura = constanteTemperatura
        self.deltaTime = 0.01
        self.centro = [self.frameSize / 2, self.frameSize / 2]
        pass

    def layout(self):
        """
        Aplica el algoritmo de Fruchtermann-Reingold para obtener (y mostrar)
        un layout
        """
        self.algoritmoFruchtermanReingold()
        pass

    def posicionesAleatorias(self):
        for vertice in self.vertices:
            self.posiciones[vertice] = np.array([np.random.random_sample() * self.frameSize, np.random.random_sample() * self.frameSize])
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

        plt.pause(self.deltaTime)
        pass

    def algoritmoFruchtermanReingold(self):
        # Seteamos posiciones iniciales aleatorias
        self.posicionesAleatorias()
        self.initializeTemperature()
        start_time = time.time()
        for count in range(self.iters):
            self.step(count)

        self.plotear()
        if self.verbose:
            print("--- Tiempo de ejecucion: %s segundos ---" % (time.time() - start_time))
        plt.show()
        pass

    def step(self, count):
        self.initializeAccumulators()
        self.computeAttractionForces()
        self.computeRepulsionForces()
        self.computeGravityForces()
        self.updatePositions()
        self.updateTemperature()
        if count % self.refresh == 0:
            self.plotear()
        self.mostrarVerbosidad(count)

        pass

    def mostrarVerbosidad(self, count):
        if self.verbose:
            print("\n Iteracion: ", count)
            print("Temperatura: ", self.temperatura)
            print("Fuerza acumulada en x: ")
            pprint.pp(self.accumx)
            print("Fuerza acumulada en y: ")
            pprint.pp(self.accumy)

    def initializeTemperature(self):
        self.temperatura = self.temperaturaInicial
        pass

    def initializeAccumulators(self):
        for vertice in self.vertices:
            self.accumx[vertice] = 0
            self.accumy[vertice] = 0
        pass

    def computeAttractionForces(self):
        for ni, nj in self.aristas:
            distance = distanciaEuclidiana(self.posiciones[ni], self.posiciones[nj])
            modfa = self.attraction(distance)
            fx = modfa * (self.abcisa(nj) - self.abcisa(ni)) / distance
            fy = modfa * (self.ordenada(nj) - self.ordenada(ni)) / distance

            self.accumx[ni] += fx
            self.accumx[nj] -= fx

            self.accumy[ni] += fy
            self.accumy[nj] -= fy
        pass

    def computeRepulsionForces(self):
        for ni in self.vertices:
            for nj in self.vertices:
                if ni != nj:
                    distance = distanciaEuclidiana(self.posiciones[ni], self.posiciones[nj])
                    modfa = self.repultion(distance)
                    fx = modfa * (self.abcisa(nj) - self.abcisa(ni)) / distance
                    fy = modfa * (self.ordenada(nj) - self.ordenada(ni)) / distance

                    self.accumx[ni] -= fx
                    self.accumx[nj] += fx

                    self.accumy[ni] -= fy
                    self.accumy[nj] += fy
        pass

    def computeGravityForces(self):
        for ni in self.vertices:
                distance = distanciaEuclidiana(self.posiciones[ni], self.centro)
                modfa = self.gravity
                fx = modfa * (self.centro[0] - self.abcisa(ni)) / distance
                fy = modfa * (self.centro[1] - self.ordenada(ni)) / distance

                self.accumx[ni] -= fx
                self.accumy[ni] -= fy
        
    def updatePositions(self):
        for node in self.vertices:
            f = np.array([self.accumx[node], self.accumy[node]])
            if modulo(f) > self.temperatura:
                f = productoPorEscalar(self.temperatura / modulo(f), f)
                self.accumx[node] = f[0]
                self.accumy[node] = f[1]
            self.posiciones[node][0] += self.accumx[node]
            self.posiciones[node][1] += self.accumy[node]

        pass

    def updateTemperature(self):
        self.temperatura = self.temperatura * self.constanteTemperatura
        pass

    def ordenada(self, v):
        return self.posiciones[v][1]

    def repultion(self, distance):
        return self.k ** 2 / distance * self.repultionConstant

    def attraction(self, distance):
        return distance ** 2 / self.k * self.attractionConstant


    def abcisa(self, v):
        return self.posiciones[v][0]


def distanciaEuclidiana(ni, nj):
    return np.linalg.norm(ni - nj)

def modulo(vector):
    return distanciaEuclidiana(vector, (0, 0))

def productoPorEscalar(escalar, vector):
    return escalar * vector


if __name__ == "__main__":
    for a, b, c in np.random.random((100, 3, 2)):
        BA = distanciaEuclidiana(b, a)

        AB = distanciaEuclidiana(a, b)
        BC = distanciaEuclidiana(b, c)
        AC = distanciaEuclidiana(a, c)

        assert(AB == BA)
        assert(AB >= 0)
        assert(AB + BC >= AC)