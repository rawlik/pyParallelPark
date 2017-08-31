#! /usr/bin/python
# -*- coding: utf-8 -*-

from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
import sys
from math import *
from random import random


old_t = 0
name = 'parkowanie'

def circ(x, y, r, ile = 10):
    """Rysowanie koła"""
    p = ceil(-log((2 * pi / ile), 10 ))
    t = []
    for phi in range(0, 2 * pow(10,p) * pi, 1):
        phi = phi / float(pow(10,p))
        t.append( (r * sin(phi) + x, r * cos(phi) + y) )
    return t

class Point:
    """Punkt który zwraca swoje położenie w układzie
    samochodu lub świata. Pamięta współrzędne w świecie!"""
    def __init__(self, car, coord = None, flag = 'world'):
        self.car = car
        if flag == 'world':
            self.SetWorld(coord)
        elif flag == 'car':
            self.SetCar(coord)
        else:
            self.coord = None

    def GetCar(self):
        x = (self.coord[0] - self.car.coord[0]) * cos(self.car.phi)
        x += (self.coord[1] - self.car.coord[1]) * sin(self.car.phi)
        y = - (self.coord[0] - self.car.coord[0]) * sin(self.car.phi)
        y += (self.coord[1] - self.car.coord[1]) * cos(self.car.phi)
        return (x, y)
    
    def GetWorld(self):
        return self.coord

    def GetOther(self, o, phi):
        x = (self.coord[0] - o[0]) * cos(phi)
        x += (self.coord[1] - o[1]) * sin(phi)
        y = - (self.coord[0] - o[0]) * sin(phi)
        y += (self.coord[1] - o[1]) * cos(phi)
        return (x, y)
    
    def SetCar(self, coord):
        x = coord[0] * cos(self.car.phi) - coord[1] * sin(self.car.phi) + self.car.coord[0]
        y = coord[0] * sin(self.car.phi) + coord[1] * cos(self.car.phi) + self.car.coord[1]
        self.coord = (x, y)
    
    def SetWorld(self, coord):
        self.coord = coord
        
    def Vertex(self):
        if self.coord != None:
            glVertex2d(self.coord[0], self.coord[1])

    def draw(self):
        glColor3d(1,0,0)
        glPointSize(5)
        glBegin(GL_POINTS)
        self.Vertex()
        glEnd()

    def __add__(self, pnt):
        return Point(self.car, (self.coord[0] + pnt[0], self.coord[1] + pnt[1]), 'world')

    def __sub__(self, pnt):
        return sqrt( (self.coord[0] - pnt.coord[0])**2 + (self.coord[1] - pnt.coord[1])**2 )

    def __getitem__(self, key):
        return self.coord[key]


class Meter:
    def __init__(self, car, coord = (0.04, -0.04)):
        self.coord = coord
        self.car = car
        self.point = Point(self.car)
        self.UpdateCoord()
        self.length = 2
        self.polygons=[]
        self.hit = False
        self.hitpoint = None
        self.read = 0

    def UpdateCoord(self):
        self.point.SetCar(self.coord)

    def GetCoord(self):
        return self.point.GetWorld()

    def GetAngle(self):
        return self.car.phi

    def GetA(self):
        """A, jesli promien miernika przebiega wg. AX+B"""
        return tan(self.GetAngle())

    def GetB(self):
        """B, jesli promien miernika przebiega wg. AX+B"""
        return self.GetCoord()[1] - self.GetA() * self.GetCoord()[0]

    def draw(self):
        """rysowanie kreski miernika,
        miernik zmienia kolor, gdy w cos trafia"""
        if self.hit:
            glColor3d(0, 1, 1)
        else:
            glColor3d(0.5, 0, 0.5)
        glBegin(GL_LINES)
        self.point.Vertex()
        if self.hit:
            self.hitpoint.Vertex()
        else:
            (self.point + (self.length * cos(self.GetAngle()), self.length * sin(self.GetAngle()))).Vertex()
        glEnd()

    def GetRead(self):
        if self.hit:
            return self.read
        else:
            return 'Inf'

    def TestIntersection(self, list):
        """testuje przeciecie z wielokatem, ktorego
        kolejne wierzcholki znajduja sie w list"""
        if list == None or len(list) == 1 or len(list) == 0:
            return False
        #punkty nad i pod prosta
        over = []
        under = []
        for i in range(len(list)):
            if list[i].GetOther(self.GetCoord(), self.GetAngle())[0] < 0:
                continue
            if self.GetA() * list[i].GetWorld()[0] - list[i].GetWorld()[1] + self.GetB() > 0:
                over.append(i)
            else:
                under.append(i)
        if len(over) == 0 or len(under) == 0:
            return (False, None, None)
        else:
            #TRAFIENIE
            #szukamy przeciec
            #lista punktów przeciec
            points = []
            for i in under:
                if over.count((i-1) % len(list)) > 0:
                    #prosta laczaca dwa punkty
                    if (list[i-1].GetWorld()[0] - list[i].GetWorld()[0]) == 0:
                        x = list[i].GetWorld()[0]
                        y = self.GetA() * x + self.GetB()
                    else:
                        a = (list[i-1].GetWorld()[1] - list[i].GetWorld()[1]) / (list[i-1].GetWorld()[0] - list[i].GetWorld()[0])
                        b = -list[i].GetWorld()[0] * a + list[i].GetWorld()[1]
                        #punkt przeciecia
                        x = (b - self.GetB()) / (self.GetA() - a)
                        y = a * x + b
                    points.append(Point(self.car, (x,y)))
                if over.count((i+1) % len(list)) > 0:
                    #prosta laczaca dwa punkty
                    if (list[(i+1) % len(list)].GetWorld()[0] - list[i].GetWorld()[0]) == 0:
                        x = list[i].GetWorld()[0]
                        y = self.GetA() * x + self.GetB()
                    else:
                        a = (list[(i+1) % len(list)].GetWorld()[1] - list[i].GetWorld()[1]) / (list[(i+1) % len(list)].GetWorld()[0] - list[i].GetWorld()[0])
                        b = -list[i].GetWorld()[0] * a + list[i].GetWorld()[1]
                        #punkt przeciecia
                        x = (b - self.GetB()) / (self.GetA() - a)
                        y = a * x + b
                    points.append(Point(self.car, (x,y)))
            dist = 10e9
            hit_point = None
            for p in points:
                dist2 = sqrt((p.GetWorld()[0] - self.GetCoord()[0])**2 + (p.GetWorld()[1] - self.GetCoord()[1])**2)
                
                if dist2 < dist:
                    dist = dist2
                    hit_point = p

            return (True, dist, hit_point)

    def Test(self):
        """testuje przeciecia ze wszystkimi wielokatami w 
        liscie polygons"""
        self.hit = False
        dist = 10e9
        for lst in self.polygons:
            ret = self.TestIntersection(lst)
            if ret[0] == True:
                self.hit = True
                if ret[1] < dist:
                    dist = ret[1]
                    self.hitpoint = ret[2]
                    self.read = ret[1]
            
        

class Car:
    """phi - kąt pod jakim samochód jedzie,
    alpha - kąt pod jakim ustawione są koła"""
    def __init__(self):
        #parametry operacyjne
        self.v = 0
        self.phi = 0
        self.alpha = 0
        self.coord = (0, 0)
        self.last_r = 0
        self.draw_additional_info = False
        self.given_distance = False

        #parametry do zmieniania
        self.width = 0.1
        self.length = 0.2
        self.speed_step = 4e-7
        self.alpha_step = 2e-3
        self.base_speed = 7e-5
        self.wheel_radius = 0.03
        self.body_width = 0.11
        self.body_length = 0.3


        #obraz samochodu
        self.chassis = []
        for i in range(6):
            self.chassis.append(Point(self))

        self.wheels = []
        for i in range(8):
            self.wheels.append(Point(self))

        self.body = []
        for i in range(4):
            self.body.append(Point(self))

        #czujnik odleglosci
        self.meter = Meter(self)

    def update_points(self):
        #dolna oś
        self.chassis[0].SetCar((-self.width / 2, 0))
        self.chassis[1].SetCar((self.width / 2, 0))
        #górna oś
        self.chassis[2].SetCar((-self.width / 2, self.length))
        self.chassis[3].SetCar((self.width / 2, self.length))
        #wał
        self.chassis[4].SetCar((0, 0))
        self.chassis[5].SetCar((0, self.length))

        #karoseria
        self.body[0].SetCar((-self.body_width / 2, -self.body_length / 2 + self.length / 2))
        self.body[1].SetCar((-self.body_width / 2, self.body_length / 2 + self.length / 2))
        self.body[2].SetCar((self.body_width / 2, self.body_length / 2 + self.length / 2))
        self.body[3].SetCar((self.body_width / 2, -self.body_length / 2 + self.length / 2))


        #lewe tylne koło
        self.wheels[0].SetCar((-self.width / 2, -self.wheel_radius))
        self.wheels[1].SetCar((-self.width / 2, self.wheel_radius))
        #prawe tylne koło
        self.wheels[2].SetCar((self.width / 2, -self.wheel_radius))
        self.wheels[3].SetCar((self.width / 2, self.wheel_radius))

        self.wheels[4].SetCar((-self.width / 2 + self.wheel_radius * sin(self.alpha), self.length - self.wheel_radius * cos(self.alpha)))
        self.wheels[5].SetCar((-self.width / 2 - self.wheel_radius * sin(self.alpha), self.length + self.wheel_radius * cos(self.alpha)))

        self.wheels[6].SetCar((self.width / 2 + self.wheel_radius * sin(self.alpha), self.length - self.wheel_radius * cos(self.alpha)))
        self.wheels[7].SetCar((self.width / 2 - self.wheel_radius * sin(self.alpha), self.length + self.wheel_radius * cos(self.alpha)))


    def move(self, t):
        #nowe polozenie po czasie t w ukladzie aktualnego polozenia
        #ruch jest po okregu o promieniu r o kat theta 

        #jeśli mamy zadaną odległość
        if self.given_distance:
            if self.distance_left < abs(self.v * t):
                #trik: cmp(x,0) zwraca znak x
                d = cmp(self.v,0) * self.distance_left
            else:
                d = self.v * t
        else:
            d = self.v * t

        new = Point(self)
        theta = 0
        if abs(self.alpha) < 0.001:
            #jedziemy prosto
            new.SetCar((0, d))
            theta = 0
            self.last_r = 0
        else:
            #skrecamy - jedziemy po okregu o promieniu r o kat theta
            r = -self.length / tan(self.alpha)
            theta = cos(self.alpha) * d / r
            new.SetCar((r - r * cos(theta), r * sin(theta)))
            self.last_r = r
        self.coord = new.GetWorld()
        self.phi -= theta
        self.update_points()
        self.meter.UpdateCoord()
        self.meter.Test()
        if self.given_distance:
            #znowu trik - znak v mowi, w ktora strone jedziemy
            self.distance_left -= cmp(self.v,0) * d
            if(self.distance_left <= 0):
                self.given_distance = False
                self.v = 0

            
    def move_distance(self, d, reverse = False):
        self.given_distance = True
        self.distance_left = d
        if reverse:
            self.v = -self.base_speed
        else:
            self.v = self.base_speed

    def inc_speed(self, t):
        self.v += self.speed_step * t

    def dec_speed(self, t):
        self.v -= self.speed_step * t

    def stop(self):
        self.v = 0

    def turn_left(self, t):
        self.alpha += self.alpha_step * t

    def turn_right(self, t):
        self.alpha -= self.alpha_step * t

    def reset(self):
        self.v = 0
        self.phi = 0
        self.alpha = 0
        self.coord = (0, 0)

    def getTurnPoint(self):
        "zwraca punkt wokół którego skręca samochód"
        if self.alpha == 0:
            return None
        else:
            return Point(self, (-self.length / tan(self.alpha), 0), 'car')

    def draw(self):
        #rysuje promien po ktorym powinien sie poruszac
        if self.draw_additional_info:
            glColor3d(0.6, 0.6, 0)
            zero = Point(self)
            zero.SetCar((0, 0))
            mid = Point(self)
            mid.SetCar((self.last_r, 0))
            top = Point(self)
            top.SetCar((0, self.length))
            glBegin(GL_LINES)
            zero.Vertex()
            mid.Vertex()
            mid.Vertex()
            top.Vertex()
            glEnd()

            t = circ(self.last_r, 0, self.last_r)
            glBegin(GL_LINE_LOOP)
            for p in t:
                pnt = Point(self, p, 'car')
                pnt.Vertex()
            glEnd()

            t = circ(self.last_r, 0, sqrt(self.last_r**2 + self.length**2))
            glBegin(GL_LINE_LOOP)
            for p in t:
                pnt = Point(self, p, 'car')
                pnt.Vertex()
            glEnd()

        glColor3d(0, 0, 1)
        glBegin(GL_LINES)
        for p in self.chassis:
            p.Vertex()
        glEnd()

        glColor3d(0, 1, 0)
        glBegin(GL_LINES)
        for p in self.wheels:
            p.Vertex()
        glEnd()

        glColor3d(1, 0, 0)
        glBegin(GL_LINE_LOOP)
        for p in self.body:
            p.Vertex()
        glEnd()
        

class AI():
    def __init__(self, car):
        self.car = car
        #w kolejce znajduja sie funkcje do wywolania
        #kazda funkcja sama sie z niej usuwa
        #gdy juz uzna to za stosowne
        #self.queue = [self.scan1, self.scan2, self.turn1, self.move1, self.turn2, self.move2, self.finish]
        self.queue = []
        self.scan_step = 0.01
        self.gap_length = 0
        self.dist_to_car = 0

    def main_loop(self):
        if self.car.given_distance:
            #jeśli kazaliśmy się ruszać to czekamy na koniec ruchu
            return
        if len(self.queue) == 0:
            return
        self.queue[0]()

    def move(self, dist = 0.1):
        self.queue.append(('move', dist))

    def setWheels(self, angle):
        self.queue.append(('turn', angle))

    def getRead(self):
        return self.car.meter.GetRead()

    def drawSafeZone(self):
        "rysuje obszar w którym samochód będzie się poruszał podczas skrętu"
        center = self.car.getTurnPoint()
        if center == None:
            return
        def eval_dist(x):
            return center - x
        dist_tab = map(eval_dist, self.car.body)
        furthest = self.car.body[ dist_tab.index(max(dist_tab)) ]

        #trik: cmp(x,0) zwraca znak x
        t = circ(center[0], center[1], center.GetCar()[0] + self.car.body_width / 2 * cmp(self.car.alpha, 0))
        glColor3d(0,1,0)
        glBegin(GL_LINE_LOOP)
        for p in t:
            pnt = Point(self, p)
            pnt.Vertex()
        glEnd()

        t = circ(center[0], center[1], furthest - center)
        glColor3d(0,1,0)
        glBegin(GL_LINE_LOOP)
        for p in t:
            pnt = Point(self, p)
            pnt.Vertex()
        glEnd()



    #kolejne polecenia
    def scan1(self):
        """dojeżdżamy do końca pierwszego auta"""
        if self.car.meter.GetRead() == 'Inf':
            self.queue.remove(self.scan1)
            return
        self.car.move_distance(self.scan_step)
         
    def scan2(self):
        """przejeżdżamy miejsce do parkowania"""
        if self.car.meter.GetRead() != 'Inf':
            self.dist_to_car = self.car.meter.GetRead()
            self.queue.remove(self.scan2)
            return
        self.car.move_distance(self.scan_step)
        self.gap_length += self.scan_step

    def turn1(self):
        self.car.alpha = -pi/4
        self.queue.remove(self.turn1)

    def move1(self):
        self.car.move_distance(self.gap_length / 2, reverse = True)
        self.queue.remove(self.move1)

    def turn2(self):
        self.car.alpha = pi/4
        self.queue.remove(self.turn2)

    def move2(self):
        self.car.move_distance(self.gap_length / 2, reverse = True)
        self.queue.remove(self.move2)

    def finish(self):
        self.car.alpha = 0
        self.queue.remove(self.finish)
        

def DrawInfo():
    glColor3d(1, 1, 1)
    write(-0.99, 0.97, "Program symulujacy ruch samochodu.    RUCH: strzalki    ZATRZYMANIE: s    RESET: spacja    DODATKOWE RYSUNKI: i   WYJSCIE: q lub esc")

    write(-0.99, -0.90, "MIER. ODLEGLOSCI   " + repr(car.meter.GetRead()))
    write(-0.99, -0.93, "PREDKOSC                 " + repr(car.v * 1e4))
    write(-0.99, -0.96, "SKR. KOL                     " + repr(car.alpha * 180 / pi))
    write(-0.99, -0.99, "PROMIEN SKRETU      " + repr(car.last_r * 10))


def write(x, y, string):
    glRasterPos2d(x, y)
    for ch in string:
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, ord(ch))

class OtherCar:
    """Klasa samochodow, ktore poprostu
    sobie stoja"""
    def __init__(self, car, coord, flag = 'horizontal'):
        self.pnts = []    
        self.pnts.append(Point(car, coord, 'world'))
        if flag == 'horizontal':
            self.pnts.append(Point(car, coord, 'world') + (0, car.body_width))
            self.pnts.append(Point(car, coord, 'world') + (car.body_length, car.body_width))
            self.pnts.append(Point(car, coord, 'world') + (car.body_length, 0))

    def draw(self):
        glColor3d(1, 1, 1)    
        glBegin(GL_LINE_LOOP)
        for p in self.pnts:
            p.Vertex()
        glEnd()


        

        
#definujemy auto i ustawiamy na poczatek
car = Car()
car.coord = (-0.25, 0.25)
car.phi = -pi/2

ai = AI(car)
draw_others = True
draw_meter = True
draw_hitpoint = True
draw_safezone = True
others = [ OtherCar(car, (-0.55, 0)), OtherCar(car, (0.25, 0)) ]

for o in others:
    car.meter.polygons.append(o.pnts)



keyStates = {'up' : False, 'down' : False, 'left': False, 'right' : False}

def main():
    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
    glutInitWindowSize(700,700)
    glutCreateWindow(name)

    glClearColor(0.,0.,0.,1.)

    glutDisplayFunc(display)
    glutIdleFunc(idle)

    glMatrixMode(GL_PROJECTION)
    gluOrtho2D(-1,1,-1,1)
    glutSwapBuffers()

    glutSpecialFunc(special_keys)
    glutSpecialUpFunc(special_up_keys)
    glutKeyboardFunc(keys)

    glutMainLoop()
    return


def keys(key, x, y):
    if key == ' ':
        car.reset()
    elif key == 'q' or key == '':
        exit()
    elif key == 'i':
        car.draw_additional_info = not car.draw_additional_info
    elif key == 's':
        car.stop()
    elif key == 't':
        ai.setWheels(90 * random() - 45)
    elif key == 'x':
        ai.move()
    elif key == 'X':
        print ai.getRead()
    elif key == 'o':
        global draw_others
        draw_others = not draw_others
    elif key == 'm':
        global draw_meter
        draw_meter = not draw_meter
    elif key == 'p':
        global draw_hitpoint
        draw_hitpoint = not draw_hitpoint
    elif key == 'z':
        global draw_safezone
        draw_safezone = not draw_safezone
    glutPostRedisplay()


def special_keys(key, x, y):
    if key == GLUT_KEY_UP:
        keyStates['up'] = True
    elif key == GLUT_KEY_DOWN:
        keyStates['down'] = True
    elif key == GLUT_KEY_LEFT:
        keyStates['left'] = True
    elif key == GLUT_KEY_RIGHT:
        keyStates['right'] = True

def special_up_keys(key, x, y):
    if key == GLUT_KEY_UP:
        keyStates['up'] = False
    elif key == GLUT_KEY_DOWN:
        keyStates['down'] = False
    elif key == GLUT_KEY_LEFT:
        keyStates['left'] = False
    elif key == GLUT_KEY_RIGHT:
        keyStates['right'] = False

def keysActions(t):
    if keyStates['up'] == True:
        car.inc_speed(t)
    if keyStates['down'] == True:
        car.dec_speed(t)
    if keyStates['left'] == True:
        car.turn_left(t)
    if keyStates['right'] == True:
        car.turn_right(t)
    glutPostRedisplay()


def idle():
    global old_t
    t = glutGet(GLUT_ELAPSED_TIME)
    passed = t - old_t
    old_t = t
    #print "PASSED: ", passed
    car.move(passed)
    ai.main_loop()
    keysActions(passed)
    glutPostRedisplay()
    #print "v =", car.v
    #print "alpha =", car.alpha * 180 / pi
    return




def display():
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    car.draw()
    DrawInfo()
    if draw_others:
        for c in others:
            c.draw()
    if draw_meter:
        car.meter.draw()
    if draw_hitpoint and car.meter.hit:
        car.meter.hitpoint.draw()
    if draw_safezone:
        ai.drawSafeZone()
    glutSwapBuffers()
    return



if __name__ == '__main__': main()
