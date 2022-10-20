#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import os
import cg_algorithms as alg
import numpy as np
from PIL import Image


if __name__ == '__main__':
    #input_file = sys.argv[1]
    #output_dir = sys.argv[2]
    input_file = 'test.txt'
    output_dir = 'res'
    os.makedirs(output_dir, exist_ok=True)

    item_dict = {}
    point_dict = {}
    pen_color = np.zeros(3, np.uint8)
    width = 0
    height = 0

    with open(input_file, 'r') as fp:
        line = fp.readline()
        while line:
            line = line.strip().split(' ')
            if line[0] == 'resetCanvas':
                width = int(line[1])
                height = int(line[2])
                item_dict = {}
                point_dict = {}
            elif line[0] == 'saveCanvas':
                save_name = line[1]
                canvas = np.zeros([height, width, 3], np.uint8)
                canvas.fill(255)
                for item_id in item_dict.keys():
                    item_type = item_dict[item_id][0]
                    p_list = item_dict[item_id][1]
                    algorithm = item_dict[item_id][2]
                    color = item_dict[item_id][3]
                    if item_type == 'line':
                        pixels = alg.draw_line(p_list, algorithm)
                    elif item_type == 'polygon':
                        pixels = alg.draw_polygon(p_list, algorithm)
                    elif item_type == 'ellipse':
                        pixels = alg.draw_ellipse(p_list)

                    point_dict[item_id] = [pixels, color]

                for pixels, color in point_dict.values():
                    for x, y in pixels:
                        canvas[height - 1 - y, x] = color
                Image.fromarray(canvas).save(os.path.join(output_dir, save_name + '.bmp'), 'bmp')

            elif line[0] == 'setColor':
                pen_color[0] = int(line[1])
                pen_color[1] = int(line[2])
                pen_color[2] = int(line[3])
            elif line[0] == 'drawLine':
                item_id = line[1]
                x0 = int(line[2])
                y0 = int(line[3])
                x1 = int(line[4])
                y1 = int(line[5])
                algorithm = line[6]
                item_dict[item_id] = ['line', [[x0, y0], [x1, y1]], algorithm, np.array(pen_color)]
            elif line[0] == 'drawPolygon':
                item_id = line[1]
                line_len = len(line)
                i = 2
                pointlist = []
                while i < line_len - 1:
                    x = int(line[i])
                    y = int(line[i + 1])
                    i += 2
                    pointlist.append([x, y])
                algorithm = line[line_len - 1]
                item_dict[item_id] = ['polygon', pointlist, algorithm, np.array(pen_color)]

            elif line[0] == 'drawEllipse':
                item_id = line[1]
                x0 = int(line[2])
                y0 = int(line[3])
                x1 = int(line[4])
                y1 = int(line[5])
                algorithm = ''
                item_dict[item_id] = ['ellipse', [[x0, y0], [x1, y1]], algorithm, np.array(pen_color)]

            elif line[0] == 'rotate':
                item_id = line[1]
                pixels = item_dict[item_id][1]
                x = int(line[2])
                y = int(line[3])
                a = int(line[4])
                res = alg.rotate(pixels, x, y, a)
                item_dict[item_id][1] = res
                
            elif line[0] == 'translate':
                item_id = line[1]
                pixels = item_dict[item_id][1]
                x = int(line[2])
                y = int(line[3])
                res = alg.translate(pixels, x, y)
                item_dict[item_id][1] = res
            ...

            line = fp.readline()

