#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

import cg_algorithms as alg
from typing import Optional
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    qApp,
    QGraphicsScene,
    QGraphicsView,
    QGraphicsItem,
    QListWidget,
    QHBoxLayout,
    QWidget,
    QStyleOptionGraphicsItem,
    QColorDialog)
from PyQt5.QtGui import QPainter, QMouseEvent, QColor
from PyQt5.QtCore import QRectF

from PyQt5 import QtCore

import math

def get_length2(a, b):
    len2 = (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1])
    return len2

# retrived from https://stackoverflow.com/questions/35176451/python-code-to-calculate-angle-between-three-point-using-their-3d-coordinates
def get_angle(a, b, c):
    ang = math.degrees(math.atan2(c[1]-b[1], c[0]-b[0]) - math.atan2(a[1]-b[1], a[0]-b[0]))
    return ang


class MyCanvas(QGraphicsView):
    """
    画布窗体类，继承自QGraphicsView，采用QGraphicsView、QGraphicsScene、QGraphicsItem的绘图框架
    """
    def __init__(self, *args):
        super().__init__(*args)
        self.main_window = None
        self.list_widget = None
        self.item_dict = {}
        self.selected_id = ''

        self.status = ''
        self.temp_algorithm = ''
        self.temp_id = ''
        self.temp_item = None

        self.me_orig = []
        self.me_dest = []

    def start_draw_line(self, algorithm, item_id):
        self.status = 'line'
        self.temp_algorithm = algorithm
        self.temp_id = item_id
        self.temp_item = None

    def start_draw_polygon(self, algorithm, item_id):
        self.status = 'polygon'
        self.temp_algorithm = algorithm
        self.temp_id = item_id
        self.temp_item = None

    def start_draw_ellipse(self, algorithm, item_id):
        self.status = 'ellipse'
        self.temp_algorithm = ''
        self.temp_id = item_id
        self.temp_item = None

    def start_draw_curve(self, algorithm, item_id):
        self.status = 'curve'
        self.temp_algorithm = algorithm
        self.temp_id = item_id
        self.temp_item = None

    def start_translate(self):
        assert self.selected_id != ''
        self.status = 'translate'
        self.temp_algorithm = ''

    def start_rotate(self):
        assert self.selected_id != ''
        self.status = 'rotate'
        self.temp_algorithm = ''
        self.me_orig = None

    def start_scale(self):
        assert self.selected_id != ''
        self.status = 'scale'
        self.temp_algorithm = ''
        self.me_orig = None

    def start_clip(self, algorithm):
        assert self.selected_id != ''
        assert self.item_dict[self.selected_id].item_type == 'line'
        self.status = 'clip'
        self.temp_algorithm = algorithm

    def finish_draw(self):
        self.temp_id = self.main_window.get_id()
        self.temp_item = None

    def clear_selection(self):
        if self.selected_id != '':
            self.item_dict[self.selected_id].selected = False
            self.selected_id = ''

    def selection_changed(self, selected):
        self.main_window.statusBar().showMessage('图元选择： %s' % selected)
        if self.selected_id != '':
            self.item_dict[self.selected_id].selected = False
            self.item_dict[self.selected_id].update()
        self.selected_id = selected
        self.item_dict[selected].selected = True
        self.item_dict[selected].update()
        self.status = ''
        self.updateScene([self.sceneRect()])

    def mousePressEvent(self, event: QMouseEvent) -> None:
        pos = self.mapToScene(event.localPos().toPoint())
        x = int(pos.x())
        y = int(pos.y())
        if self.status == 'line':
            self.temp_item = MyItem(self.temp_id, self.status, [[x, y], [x, y]], self.temp_algorithm)
            self.scene().addItem(self.temp_item)
        elif self.status == 'polygon' or self.status == 'curve':
            if self.temp_item is None:
                self.temp_item = MyItem(self.temp_id, self.status, [], self.temp_algorithm)
                self.scene().addItem(self.temp_item)
            if event.button() == QtCore.Qt.LeftButton:
                self.temp_item.p_list.append([x, y])
            elif event.button() == QtCore.Qt.RightButton:
                self.item_dict[self.temp_id] = self.temp_item
                self.list_widget.addItem(self.temp_id)
                self.finish_draw()

        elif self.status == 'ellipse':
            self.temp_item = MyItem(self.temp_id, self.status, [[x, y], [x, y]], self.temp_algorithm)
            self.scene().addItem(self.temp_item)

        elif self.status == 'translate':
            if self.selected_id != '':
                self.temp_item = self.item_dict[self.selected_id]
                self.me_orig = [x, y]
        elif self.status == 'rotate':
            if self.selected_id != '':
                self.temp_item = self.item_dict[self.selected_id]
                if event.button() == QtCore.Qt.RightButton:
                    self.me_orig = [x, y]
                elif event.button() == QtCore.Qt.LeftButton:
                    self.me_dest = [x, y]
        elif self.status == 'scale':
            if self.selected_id != '':
                self.temp_item = self.item_dict[self.selected_id]
                if event.button() == QtCore.Qt.RightButton:
                    self.me_orig = [x, y]
                elif event.button() == QtCore.Qt.LeftButton:
                    self.me_dest = [x, y]
        elif self.status == 'clip':
            if self.selected_id != '':
                self.temp_item = self.item_dict[self.selected_id]
                self.me_orig = [x, y]
        elif self.status == 'selecting':
            pass
        self.updateScene([self.sceneRect()])
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:
        pos = self.mapToScene(event.localPos().toPoint())
        x = int(pos.x())
        y = int(pos.y())
        if self.status == 'line':
            self.temp_item.p_list[1] = [x, y]
        elif self.status == 'polygon' or self.status == 'curve':
            self.temp_item.p_list[-1] = [x, y]
        elif self.status == 'ellipse':
            self.temp_item.p_list[1] = [x, y]
        elif self.status == 'translate':
            self.temp_item.p_list = alg.translate(self.temp_item.p_list, x - self.me_orig[0], y - self.me_orig[1])
            self.me_orig = [x, y]
        elif self.status == 'rotate':
            if self.me_orig is not None:
                r = get_angle(self.me_dest, self.me_orig, [x, y])
                if abs(r) > 1:
                    self.me_dest = [x, y]
                    print(r)
                    self.temp_item.p_list = alg.rotate(self.temp_item.p_list, self.me_orig[0], self.me_orig[1], r)
        elif self.status == 'scale':
            if self.me_orig is not None:
                s = get_length2(self.me_orig, [x, y]) / get_length2(self.me_orig, self.me_dest)
                if s > 1.1 or 1 / s > 1.1:
                    self.me_dest = [x, y]
                    print(s)
                    self.temp_item.p_list = alg.scale(self.temp_item.p_list, self.me_orig[0], self.me_orig[1], s)
        elif self.status == 'clip':
            # TODO: draw bounding rect
            self.me_dest = [x, y]
            x0 = min(self.me_orig[0], x)
            x1 = max(self.me_orig[0], x)
            y0 = min(self.me_orig[1], y)
            y1 = max(self.me_orig[1], y)
            self.temp_item.p_list = alg.clip(self.temp_item.p_list, x0, y0, x1, y1, self.temp_algorithm)

        self.updateScene([self.sceneRect()])
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        if self.status == 'line':
            self.item_dict[self.temp_id] = self.temp_item
            self.list_widget.addItem(self.temp_id)
            self.finish_draw()
        elif self.status == 'polygon' or self.status == 'curve':
            pass
        elif self.status == 'ellipse':
            self.item_dict[self.temp_id] = self.temp_item
            self.list_widget.addItem(self.temp_id)
            self.finish_draw()
        elif self.status == 'translate':
            pass
        elif self.status == 'rotate':
            pass
        elif self.status == 'scale':
            pass
        elif self.status == 'clip':
            x0 = min(self.me_orig[0], self.me_dest[0])
            x1 = max(self.me_orig[0], self.me_dest[0])
            y0 = min(self.me_orig[1], self.me_dest[1])
            y1 = max(self.me_orig[1], self.me_dest[1])
            self.temp_item.p_list = alg.clip(self.temp_item.p_list, x0, y0, x1, y1, self.temp_algorithm)

        super().mouseReleaseEvent(event)


class MyItem(QGraphicsItem):
    """
    自定义图元类，继承自QGraphicsItem
    """
    def __init__(self, item_id: str, item_type: str, p_list: list, algorithm: str = '', parent: QGraphicsItem = None):
        """

        :param item_id: 图元ID
        :param item_type: 图元类型，'line'、'polygon'、'ellipse'、'curve'等
        :param p_list: 图元参数
        :param algorithm: 绘制算法，'DDA'、'Bresenham'、'Bezier'、'B-spline'等
        :param parent:
        """
        super().__init__(parent)
        self.id = item_id           # 图元ID
        self.item_type = item_type  # 图元类型，'line'、'polygon'、'ellipse'、'curve'等
        self.p_list = p_list        # 图元参数
        self.algorithm = algorithm  # 绘制算法，'DDA'、'Bresenham'、'Bezier'、'B-spline'等
        self.selected = False

    def paint(self, painter: QPainter, option: QStyleOptionGraphicsItem, widget: Optional[QWidget] = ...) -> None:
        item_pixels = []
        if self.item_type == 'line':
            item_pixels = alg.draw_line(self.p_list, self.algorithm)
        elif self.item_type == 'polygon':
            item_pixels = alg.draw_polygon(self.p_list, self.algorithm)
        elif self.item_type == 'ellipse':
            item_pixels = alg.draw_ellipse(self.p_list)
        elif self.item_type == 'curve':
            item_pixels = alg.draw_curve(self.p_list, self.algorithm)

        for p in item_pixels:
            painter.drawPoint(*p)
        if self.selected:
            painter.setPen(QColor(255, 0, 0))
            painter.drawRect(self.boundingRect())

    def boundingRect(self) -> QRectF:
        if self.item_type == 'line':
            x0, y0 = self.p_list[0]
            x1, y1 = self.p_list[1]
            x = min(x0, x1)
            y = min(y0, y1)
            w = max(x0, x1) - x
            h = max(y0, y1) - y
            return QRectF(x - 1, y - 1, w + 2, h + 2)
        elif self.item_type == 'polygon':
            x = min([p[0] for p in self.p_list])
            y = min([p[1] for p in self.p_list])
            w = max([p[0] for p in self.p_list]) - x
            h = max([p[1] for p in self.p_list]) - y
            return QRectF(x - 1, y - 1, w + 2, h + 2)
        elif self.item_type == 'ellipse':
            x0, y0 = self.p_list[0]
            x1, y1 = self.p_list[1]
            assert [x0, y0, x1, y1] == [min(x0, x1), max(y0, y1), max(x0, x1), min(y0, y1)]
            x = min(x0, x1)
            y = min(y0, y1)
            w = max(x0, x1) - x
            h = max(y0, y1) - y
            return QRectF(x - 1, y - 1, w + 2, h + 2)
        elif self.item_type == 'curve':
            x = min([p[0] for p in self.p_list])
            y = min([p[1] for p in self.p_list])
            w = max([p[0] for p in self.p_list]) - x
            h = max([p[1] for p in self.p_list]) - y
            return QRectF(x - 1, y - 1, w + 2, h + 2)


class MainWindow(QMainWindow):
    """
    主窗口类
    """
    def __init__(self):
        super().__init__()
        self.item_cnt = 0

        # 使用QListWidget来记录已有的图元，并用于选择图元。注：这是图元选择的简单实现方法，更好的实现是在画布中直接用鼠标选择图元
        self.list_widget = QListWidget(self)
        self.list_widget.setMinimumWidth(200)

        # 使用QGraphicsView作为画布
        self.scene = QGraphicsScene(self)
        self.scene.setSceneRect(0, 0, 600, 600)
        self.canvas_widget = MyCanvas(self.scene, self)
        self.canvas_widget.setFixedSize(600, 600)
        self.canvas_widget.main_window = self
        self.canvas_widget.list_widget = self.list_widget

        # 设置菜单栏
        menubar = self.menuBar()
        file_menu = menubar.addMenu('文件')
        set_pen_act = file_menu.addAction('设置画笔')
        reset_canvas_act = file_menu.addAction('重置画布')
        exit_act = file_menu.addAction('退出')
        draw_menu = menubar.addMenu('绘制')
        line_menu = draw_menu.addMenu('线段')
        line_naive_act = line_menu.addAction('Naive')
        line_dda_act = line_menu.addAction('DDA')
        line_bresenham_act = line_menu.addAction('Bresenham')
        polygon_menu = draw_menu.addMenu('多边形')
        polygon_dda_act = polygon_menu.addAction('DDA')
        polygon_bresenham_act = polygon_menu.addAction('Bresenham')
        ellipse_act = draw_menu.addAction('椭圆')
        curve_menu = draw_menu.addMenu('曲线')
        curve_bezier_act = curve_menu.addAction('Bezier')
        curve_b_spline_act = curve_menu.addAction('B-spline')
        edit_menu = menubar.addMenu('编辑')
        translate_act = edit_menu.addAction('平移')
        rotate_act = edit_menu.addAction('旋转')
        scale_act = edit_menu.addAction('缩放')
        clip_menu = edit_menu.addMenu('裁剪')
        clip_cohen_sutherland_act = clip_menu.addAction('Cohen-Sutherland')
        clip_liang_barsky_act = clip_menu.addAction('Liang-Barsky')

        # 连接信号和槽函数
        exit_act.triggered.connect(qApp.quit)

        line_naive_act.triggered.connect(self.line_naive_action)
        line_dda_act.triggered.connect(self.line_dda_action)
        line_bresenham_act.triggered.connect(self.line_bresenham_action)

        polygon_dda_act.triggered.connect(self.polygon_dda_action)
        polygon_bresenham_act.triggered.connect(self.polygon_bresenham_action)

        ellipse_act.triggered.connect(self.ellipse_action)

        curve_bezier_act.triggered.connect(self.curve_bezier_action)
        curve_b_spline_act.triggered.connect(self.curve_b_spline_action)

        translate_act.triggered.connect(self.translate_action)
        rotate_act.triggered.connect(self.rotate_action)
        scale_act.triggered.connect(self.scale_action)

        clip_cohen_sutherland_act.triggered.connect(self.clip_cohen_sutherland_action)
        clip_liang_barsky_act.triggered.connect(self.clip_liang_barsky_action)

        self.list_widget.currentTextChanged.connect(self.canvas_widget.selection_changed)

# TODO
#         menubar.triggered.connect = self.menuBar()
#         file_menu = menubar.addMenu('文件')
#         set_pen_act = file_menu.addAction('设置画笔')
#         reset_canvas_act = file_menu.addAction('重置画布')
#         exit_act = file_menu.addAction('退出')
#         draw_menu = menubar.addMenu('绘制')
#         line_menu = draw_menu.addMenu('线段')
#         line_naive_act = line_menu.addAction('Naive')
#         line_dda_act = line_menu.addAction('DDA')
#         line_bresenham_act = line_menu.addAction('Bresenham')
#         polygon_menu = draw_menu.addMenu('多边形')
#         polygon_dda_act = polygon_menu.addAction('DDA')
#         polygon_bresenham_act = polygon_menu.addAction('Bresenham')
#         ellipse_act = draw_menu.addAction('椭圆')
#         curve_menu = draw_menu.addMenu('曲线')
#         curve_bezier_act = curve_menu.addAction('Bezier')
#         curve_b_spline_act = curve_menu.addAction('B-spline')
#         edit_menu = menubar.addMenu('编辑')
#         translate_act = edit_menu.addAction('平移')
#         rotate_act = edit_menu.addAction('旋转')
#         scale_act = edit_menu.addAction('缩放')
#         clip_menu = edit_menu.addMenu('裁剪')
#         clip_cohen_sutherland_act = clip_menu.addAction('Cohen-Sutherland')
#         clip_liang_barsky_act = clip_menu.addAction('Liang-Barsky')

        # 设置主窗口的布局
        self.hbox_layout = QHBoxLayout()
        self.hbox_layout.addWidget(self.canvas_widget)
        self.hbox_layout.addWidget(self.list_widget, stretch=1)
        self.central_widget = QWidget()
        self.central_widget.setLayout(self.hbox_layout)
        self.setCentralWidget(self.central_widget)
        self.statusBar().showMessage('空闲')
        self.resize(600, 600)
        self.setWindowTitle('CG Demo')

    def get_id(self):
        _id = str(self.item_cnt)
        self.item_cnt += 1
        return _id

    def set_pen_action(self):
        color = QColorDialog.getColor()
        if color.isValid():
            self.canvas_widget.temp_color = color

    def line_naive_action(self):
        self.canvas_widget.start_draw_line('Naive', self.get_id())
        self.statusBar().showMessage('Naive算法绘制线段')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection()

    def line_dda_action(self):
        self.canvas_widget.start_draw_line('DDA', self.get_id())
        self.statusBar().showMessage('DDA算法绘制线段')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection()

    def line_bresenham_action(self):
        self.canvas_widget.start_draw_line('Bresenham', self.get_id())
        self.statusBar().showMessage('Bresenham算法绘制线段')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection()

    def polygon_dda_action(self):
        self.canvas_widget.start_draw_polygon('DDA', self.get_id())
        self.statusBar().showMessage('DDA算法绘制多边形')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection()

    def polygon_bresenham_action(self):
        self.canvas_widget.start_draw_polygon('Bresenham', self.get_id())
        self.statusBar().showMessage('Bresenham算法绘制多边形')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection()

    def ellipse_action(self):
        self.canvas_widget.start_draw_ellipse('', self.get_id())
        self.statusBar().showMessage('中点圆算法绘制椭圆')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection()

    def curve_bezier_action(self):
        self.canvas_widget.start_draw_curve('Bezier', self.get_id())
        self.statusBar().showMessage('Bezier算法绘制曲线')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection()

    def curve_b_spline_action(self):
        self.canvas_widget.start_draw_curve('B-spline', self.get_id())
        self.statusBar().showMessage('B-spline算法绘制曲线')
        self.list_widget.clearSelection()
        self.canvas_widget.clear_selection()

    def translate_action(self):
        self.canvas_widget.start_translate()
        self.statusBar().showMessage('平移')

    def rotate_action(self):
        self.canvas_widget.start_rotate()
        self.statusBar().showMessage('旋转')

    def scale_action(self):
        self.canvas_widget.start_scale()
        self.statusBar().showMessage('缩放')

    def clip_cohen_sutherland_action(self):
        self.canvas_widget.start_clip('Cohen-Sutherland')
        self.statusBar().showMessage('Cohen-Sutherland算法裁剪线段')

    def clip_liang_barsky_action(self):
        self.canvas_widget.start_clip('Liang-Barsky')
        self.statusBar().showMessage('Liang-Barsky算法裁剪线段')


if __name__ == '__main__':
    app = QApplication(sys.argv)
    mw = MainWindow()
    mw.show()
    sys.exit(app.exec_())
