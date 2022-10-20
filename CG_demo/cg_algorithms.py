#!/usr/bin/env python
# -*- coding:utf-8 -*-

# 本文件只允许依赖math库
import math


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
        dx = x1 - x0
        dy = y1 - y0
        abs_dx = abs(dx)
        abs_dy = abs(dy)
        if abs_dx >= abs_dy:
            step = abs_dx
        else:
            step = abs_dy
        increx = float(dx)/step
        increy = float(dy)/step
        cur_x = float(x0)
        cur_y = float(y0)
        while (cur_x != float(x1)) | (cur_y != float(y1)):
            result.append((int(cur_x), int(cur_y)))
            cur_x += increx
            cur_y += increy
        result.append((x1, y1))
    elif algorithm == 'Bresenham':
        dx = abs(x1-x0)
        if x0 < x1:
            sx = 1
        else:
            sx = -1
        dy = abs(y1-y0)
        if y0 < y1:
            sy = 1
        else:
            sy = -1
         
        t = 0
        if dx > dy:
            tmp = dx
            dx = dy
            dy = tmp
            t = 1
        err = 2 * dy - dx
        while (x0 != x1) | (y0 != y1):
            result.append((x0, y0))
            e2 = err
            if e2 >= 0:
                if t == 1:
                    x0 += sx
                else:
                    y0 += sy
                err = err - 2 * dx
            else:
                if t == 1:
                    y0 += sy
                else:
                    x0 += sx
                err = err + 2 * dy

        result.append((x1, y1))
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


def draw_ellipse(p_list):
    """绘制椭圆（采用中点圆生成算法）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 椭圆的矩形包围框左上角和右下角顶点坐标
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    x0 = int(p_list[0][0])
    y0 = int(p_list[0][1])
    x1 = int(p_list[1][0])
    y1 = int(p_list[1][1])
    if abs(x1 - x0) % 2 == 1:
        a = (abs(x1 - x0) + 1) / 2
    else:
        a = abs(x1 - x0) / 2
    if abs(y1 - y0) % 2 == 1:
        b = (abs(y0 - y1) + 1) / 2
    else:
        b = abs(y0 - y1) / 2
    cur_x = 0
    cur_y = int(b)
    pl = pow(b, 2) - pow(a, 2) * b + pow(a, 2) / 4
    judge = pow(b, 2)*cur_x - pow(a, 2)*cur_y
    res = []
    while judge < 0:
        res.append((int(cur_x), int(cur_y)))
        if pl < 0:
            next_x = cur_x + 1
            next_y = cur_y
            pl = pl + 2 * pow(b, 2) * next_x + pow(b, 2)
        else:
            next_x = cur_x + 1
            next_y = cur_y - 1
            pl = pl + 2 * pow(b, 2) * next_x - 2 * pow(a, 2) * next_y + pow(b, 2)
        cur_x = next_x
        cur_y = next_y
        judge = pow(b, 2)*cur_x - pow(a, 2)*cur_y

    pl = pow(b, 2) * (cur_x + 1/2) * (cur_x + 1/2) + pow(a, 2) * (cur_y - 1) * (cur_y - 1) - pow(b, 2) * pow(a, 2)

    while cur_y > 0:
        res.append((int(cur_x), int(cur_y)))
        if pl >= 0:
            next_x = cur_x
            next_y = cur_y - 1
            pl = pl - 2 * pow(a, 2) * next_y + pow(a, 2)
        else:
            next_x = cur_x + 1
            next_y = cur_y - 1
            pl = pl + 2 * pow(b, 2) * next_x - 2 * pow(a, 2) * next_y + pow(a, 2)
        cur_x = next_x
        cur_y = next_y
    res.append((cur_x, cur_y))
    final_res = []
    for point in res:
        x = int(point[0])
        y = int(point[1])
        final_res.append((x, y))
        final_res.append((x, -y))
        final_res.append((-x, y))
        final_res.append((-x, -y))
    return translate(final_res, int((x0 + x1)/2), int((y0 + y1)/2))


def draw_curve(p_list, algorithm):
    """绘制曲线

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 曲线的控制点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'Bezier'和'B-spline'（三次均匀B样条曲线，曲线不必经过首末控制点）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    pass


def translate(p_list, dx, dy):
    """平移变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param dx: (int) 水平方向平移量
    :param dy: (int) 垂直方向平移量
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    res = []
    for point in p_list:
        x = int(point[0]) + dx
        y = int(point[1]) + dy
        res.append((x, y))
    return res


def rotate(p_list, x, y, r):
    """旋转变换（除椭圆外）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 旋转中心x坐标
    :param y: (int) 旋转中心y坐标
    :param r: (int) 顺时针旋转角度（°）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    cos = math.cos(r/180 * math.pi)
    sin = math.sin(r/180 * math.pi)
    res = []
    for point in p_list:
        cur_x = int(point[0])
        cur_y = int(point[1])
        next_x = x + (cur_x - x) * cos - (cur_y - y) * sin
        next_y = y + (cur_x - x) * sin + (cur_y - y) * cos
        res.append((int(next_x), int(next_y)))
    return res


def scale(p_list, x, y, s):
    """缩放变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 缩放中心x坐标
    :param y: (int) 缩放中心y坐标
    :param s: (float) 缩放倍数
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    pass


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
    pass
