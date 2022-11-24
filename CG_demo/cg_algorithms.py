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
    if p_list == []:
        return []
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
    res = []
    for point in p_list:
        dx = int(point[0]) - x
        dy = int(point[1]) - y
        newx = int(dx * s + x)
        newy = int(dy * s + y)
        res.append([newx, newy])
    return res


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
    x0 = int(p_list[0][0])
    y0 = int(p_list[0][1])
    x1 = int(p_list[1][0])
    y1 = int(p_list[1][1])
    if algorithm == 'Cohen-Sutherland':
        j1 = 0
        j2 = 0
        if x0 - x_min > 0:
            j1 += 1
        if x_max - x0 > 0:
            j1 += 2
        if y0 - y_min > 0:
            j1 += 4
        if y_max - y0 > 0:
            j1 += 8

        if x1 - x_min > 0:
            j2 += 1
        if x_max - x1 > 0:
            j2 += 2
        if y1 - y_min > 0:
            j2 += 4
        if y_max - y1 > 0:
            j2 += 8

        if (j1 == 0) & (j2 == 0):
            return p_list
        elif (j1 & j2) != 0:
            return []
        if x1 == x0:
            res_y0 = 0
            res_y1 = 0
            if y0 > y1:
                max_y = y0
                min_y = y1
            else:
                max_y = y1
                min_y = y0

            if max_y > y_max:
                res_y0 = y_max
            else:
                res_y0 = max_y
            if min_y < y_min:
                res_y1 = y_min
            else:
                res_y1 = min_y
            return [[x0, int(res_y0)], [x1, int(res_y1)]]
        elif y1 == y0:
            res_x0 = 0
            res_x1 = 0
            if x0 > x1:
                max_x = x0
                min_x = x1
            else:
                max_x = x1
                min_x = x0
            if max_x > x_max:
                res_x0 = x_max
            else:
                res_x0 = max_x
            if min_x < x_min:
                res_x1 = x_min
            else:
                res_x1 = min_x
            return [[int(res_x0), y0], [int(res_x1), y1]]
        else:
            if x0 > x1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            k = (y1 - y0) / (x1 - x0)
            m = (x1 - x0) / (y1 - y0)
            y_left = k * (x_min - x0) + y0
            y_right = k * (x_max - x0) + y0
            x_up = m * (y_max - y0) + x0
            x_down = m * (y_min - y0) + x0
            judge = False
            if (y_left <= y_max) & (y_left >= y_min):
                if (x0 <= x_min) & (x1 >= x_min):
                    x0 = x_min
                    y0 = y_left
                    judge = True
            if (y_right <= y_max) & (y_left >= y_min):
                if (x0 <= x_max) & (x1 >= x_max):
                    x1 = x_max
                    y1 = y_right
                    judge = True

            if y0 > y1:
               x0, y0, x1, y1 = x1, y1, x0, y0
            if (x_up <= x_max) & (x_up >= x_min):
                if (y0 <= y_max) & (y1 >= y_max):
                    x1 = x_up
                    y1 = y_max
                    judge = True
            if (x_down <= x_max) & (x_down >= x_min):
                if (y0 <= y_min) & (y1 >= y_min):
                    x0 = x_down
                    y0 = x_min
                    judge = True

            if judge:
                return [[int(x0), int(y0)], [int(x1), int(y1)]]
    elif algorithm == 'Liang-Barsky':
        if x0 > x1:
            x0, y0, x1, y1 = x1, y1, x0, y0

        dx = x1 - x0
        dy = y1 - y0
        p = [x0-x1, x1-x0, y0-y1, y1-y0]
        q = [x0-x_min, x_max-x0, y0-y_min, y_max-y0]
        u0, u1 = 0, 1
        for i in range(4):
            if p[i] < 0:
                u0 = max(u0, q[i]/p[i])
            elif p[i] > 0:
                u1 = min(u1, q[i]/p[i])
            elif (p[i] == 0) and (q[i] < 0):
                return []
            if u0 > u1:
                return []
        res_x0 = 0
        res_y0 = 0
        res_x1 = 0
        res_y1 = 0
        if u0 > 0:
            res_x0 = int(x0 + u0*(x1-x0) + 0.5)
            res_y0 = int(y0 + u0*(y1-y0) + 0.5)
        if u1 < 1:
            res_x1 = int(x0 + u1*(x1-x0) + 0.5)
            res_y1 = int(y0 + u1*(y1-y0) + 0.5)
        result = [[res_x0, res_y0], [res_x1, res_y1]]

        return result

