#!/usr/bin/env python
# -*- coding:utf-8 -*-

# 本文件只允许依赖math库
import math


def reflex_handler(x, y, reflex):
    if reflex:
        return y, x
    else:
        return x, y


def bspline_handler(u, p_list):
    ret = (-u ** 3 + 3 * u ** 2 - 3 * u + 1) * p_list[0] + \
          (3 * u ** 3 - 6 * u ** 2 + 4) * p_list[1] + \
          (-3 * u ** 3 + 3 * u ** 2 + 3 * u + 1) * p_list[2] + \
          (u ** 3) * p_list[3]
    return ret / 6


def cohen_sutherland_encoder(x, y, x_min, y_min, x_max, y_max):
    return ((y > y_max) << 3) + ((y < y_min) << 2) + ((x > x_max) << 1) + ((x < x_min) << 0)


def cohen_sutherland_handler(x0, y0, x1, y1, x):
    return [x, y0 + (y1 - y0) / (x1 - x0) * (x - x0)]


def draw_line(p_list, algorithm):
    """绘制线段

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'，此处的'Naive'仅作为示例，测试时不会出现
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    result = []
    if algorithm == 'Naive':
        if x0 == x1:
            for y in range(y0, y1 + 1):
                result.append((x0, y))
        else:
            if x0 > x1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            k = (y1 - y0) / (x1 - x0)
            for x in range(x0, x1 + 1):
                result.append((x, int(y0 + k * (x - x0))))
    elif algorithm == 'DDA':
        # wikipedia
        dx = x1 - x0
        dy = y1 - y0
        if abs(dx) >= abs(dy):
            step = abs(dx)
        else:
            step = abs(dy)
        if step == 0:
            return result
        dx = dx / step
        dy = dy / step
        if abs(dx) == 1:
            dx = int(dx)
        elif abs(dy) == 1:
            dy = int(dy)
        x = x0
        y = y0
        for i in range(0, step + 1):
            result.append((int(x), int(y)))
            x = x + dx
            y = y + dy

    elif algorithm == 'Bresenham':
        if x0 == x1:
            for y in range(min(y0, y1), max(y0, y1) + 1):
                result.append((x0, y))
        elif y0 == y1:
            for x in range(min(x0, x1), max(x0, x1) + 1):
                result.append((x, y0))
        else:
            # 扩展大于 1 斜率
            if abs(y1 - y0) > abs(x1 - x0):
                steep = True
                x0, y0, x1, y1 = y0, x0, y1, x1
            else:
                steep = False
            # 扩展反向
            if x0 > x1:
                x0, y0, x1, y1 = x1, y1, x0, y0

            dx = x1 - x0
            dy = y1 - y0
            yi = 1
            # 扩展负斜率
            if dy < 0:
                yi = -1
                dy = -dy
            d = 2 * dy - dx
            y = y0

            for x in range(x0, x1 + 1):
                result.append(reflex_handler(int(x), int(y), steep))
                if d > 0:
                    y = y + yi
                    d = d + 2 * (dy - dx)
                else:
                    d = d + 2 * dy
    return result


def draw_polygon(p_list, algorithm):
    """绘制多边形

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 多边形的顶点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    result = []
    for i in range(len(p_list)):
        line = draw_line([p_list[i - 1], p_list[i]], algorithm)
        result += line
    return result


def draw_ellipse_half(xo, yo, rx, ry, reflex):
    result = []
    x = 0
    y = ry
    p1 = ry ** 2 - rx ** 2 * ry + 0.25 * rx ** 2
    while ry * ry * x <= rx * rx * y:
        # 四个端点会重复
        ax, ay = reflex_handler(int(x), int(y), reflex)
        result.append((int(xo + ax), int(yo + ay)))
        result.append((int(xo + ax), int(yo - ay)))
        result.append((int(xo - ax), int(yo + ay)))
        result.append((int(xo - ax), int(yo - ay)))
        if p1 < 0:
            p1 = p1 + 2 * ry ** 2 * x + 3 * ry ** 2
        else:
            p1 = p1 + 2 * ry ** 2 * x - 2 * rx ** 2 * y + 2 * rx ** 2 + 3 * ry ** 2
            y = y - 1
        x = x + 1
    return result


def draw_ellipse(p_list):
    """绘制椭圆（采用中点圆生成算法）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 椭圆的矩形包围框左上角和右下角顶点坐标
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    result = []
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    if (x1 - x0) % 2:
        x1 = x1 - 1
    if (y1 - y0) % 2:
        y1 = y1 + 1

    xo, yo = (x0 + x1) / 2, (y0 + y1) / 2
    rx, ry = (x1 - x0) / 2, (y0 - y1) / 2

    result += draw_ellipse_half(xo, yo, rx, ry, False)
    result += draw_ellipse_half(xo, yo, ry, rx, True)

    return result


def draw_curve(p_list, algorithm):
    """绘制曲线

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 曲线的控制点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'Bezier'和'B-spline'（三次均匀B样条曲线，曲线不必经过首末控制点）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    result = []
    if algorithm == 'Bezier':
        du = 0.001
        result.append(p_list[0])
        lx, ly = p_list[0]
        u = du
        while u < 1:
            inter = p_list.copy()
            for i in range(len(p_list) - 1):
                for j in range(len(inter) - 1):
                    inter[j] = (1 - u) * inter[j][0] + u * inter[j + 1][0], (1 - u) * inter[j][1] + u * inter[j + 1][1]
            if [lx, ly] != [int(inter[0][0]), int(inter[0][1])]:
                lx, ly = int(inter[0][0]), int(inter[0][1])
                result.append([lx, ly])
            u = u + du
        result.append(p_list[-1])
    elif algorithm == 'B-spline':
        du = 0.001
        u = 0
        while u <= 1:
            for i in range(len(p_list) - 3):
                x = bspline_handler(u, [point[0] for point in p_list[i:i + 4]])
                y = bspline_handler(u, [point[1] for point in p_list[i:i + 4]])
                result.append([int(x), int(y)])
            u = u + du
    return result


def translate(p_list, dx, dy):
    """平移变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param dx: (int) 水平方向平移量
    :param dy: (int) 垂直方向平移量
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    return [[p[0] + dx, p[1] + dy] for p in p_list]


def rotate(p_list, x, y, r):
    """旋转变换（除椭圆外）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 旋转中心x坐标
    :param y: (int) 旋转中心y坐标
    :param r: (int) 顺时针旋转角度（°）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    theta = r * math.pi / 180
    return [[int((p[0] - x) * math.cos(theta) - (p[1] - y) * math.sin(theta) + x),
             int((p[0] - x) * math.sin(theta) + (p[1] - y) * math.cos(theta) + y)] for p in p_list]


def scale(p_list, x, y, s):
    """缩放变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 缩放中心x坐标
    :param y: (int) 缩放中心y坐标
    :param s: (float) 缩放倍数
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    return [[int((p[0] - x) * s + 0.5) + x, int((p[1] - y) * s + 0.5) + y] for p in p_list]


def clip(p_list, x_min, y_min, x_max, y_max, algorithm):
    """线段裁剪

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param x_min: 裁剪窗口左上角x坐标
    :param y_min: 裁剪窗口左上角y坐标
    :param x_max: 裁剪窗口右下角x坐标
    :param y_max: 裁剪窗口右下角y坐标
    :param algorithm: (string) 使用的裁剪算法，包括'Cohen-Sutherland'和'Liang-Barsky'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1]]) 裁剪后线段的起点和终点坐标
    """
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    if algorithm == 'Cohen-Sutherland':
        while True:
            code0 = cohen_sutherland_encoder(x0, y0, x_min, y_min, x_max, y_max)
            code1 = cohen_sutherland_encoder(x1, y1, x_min, y_min, x_max, y_max)

            if code0 & code1 != 0:
                return []
            elif code0 | code1 == 0:
                return [[int(x0), int(y0)], [int(x1), int(y1)]]
            else:
                if code0 == 0:
                    x0, y0, code0, x1, y1, code1 = x1, y1, code1, x0, y0, code0,
                if (code0 >> 3) & 1:
                    y0, x0 = cohen_sutherland_handler(y0, x0, y1, x1, y_max)
                if (code0 >> 2) & 1:
                    y0, x0 = cohen_sutherland_handler(y0, x0, y1, x1, y_min)
                if (code0 >> 1) & 1:
                    x0, y0 = cohen_sutherland_handler(x0, y0, x1, y1, x_max)
                if (code0 >> 0) & 1:
                    x0, y0 = cohen_sutherland_handler(x0, y0, x1, y1, x_min)

    elif algorithm == 'Liang-Barsky':
        dx = x1 - x0
        dy = y1 - y0
        p = [-dx, dx, -dy, dy]
        q = [x0 - x_min, x_max - x0, y0 - y_min, y_max - y0]
        u1, u2 = 0.0, 1.0
        for i in range(4):
            if p[i] == 0 and q[i] < 0:
                return []
            elif p[i] < 0: # outside to inside
                u1 = max(u1, q[i] / p[i])
            elif p[i] > 0: # inside to outside
                u2 = min(u2, q[i] / p[i])
            if u1 > u2:
                return []
        return [[int(x0 + u1 * dx), int(y0 + u1 * dy)],
                [int(x0 + u2 * dx), int(y0 + u2 * dy)]]
