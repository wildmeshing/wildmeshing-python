import json
import math
import wildmeshing.parse_svg.svgpathtools as svg
import numpy as np

from wildmeshing import triangulate_data

def triangulate_svg(svg_path,
        stop_quality = -1,
        max_its = 80,
        stage = 1,
        epsilon = -1,
        feature_epsilon = 1e-3,
        target_edge_len = -1,
        edge_length_r = 1./20.,
        flat_feature_angle = 10.,
        cut_outside = False,
        skip_eps = False,
        hole_file = "",
        mute_log = False):
    vertices, lines, json = convert_svg(svg_path)

    V_out, F_out, nodes, F_nodes = triangulate_data(vertices, lines, json,
    stop_quality,  max_its,  stage,
         epsilon,  feature_epsilon,
         target_edge_len,  edge_length_r,
         flat_feature_angle,
         cut_outside,
         skip_eps,
         hole_file,
         mute_log
    )

    return V_out, F_out, nodes, F_nodes






class RationalBezier(object):
    def __init__(self, arc, tstart=0, tend=1):
        self.start = arc.point(tstart)
        self.end = arc.point(tend)

        s_tangentc = arc.unit_tangent(tstart)
        e_tangentc = arc.unit_tangent(tend)


        s_tangent = (arc.tf.dot(  np.array([[s_tangentc.real], [s_tangentc.imag], [0.0]])))
        e_tangent = (arc.tf.dot(  np.array([[e_tangentc.real], [e_tangentc.imag], [0.0]])))

        self.weights = np.array([1., math.cos((tend-tstart)*arc.delta/180.*math.pi/2.), 1.])

        sx = self.start.real
        sy = self.start.imag

        ex = self.end.real
        ey = self.end.imag

        stx = s_tangent[0]
        sty = s_tangent[1]
        stn = math.sqrt(stx**2+sty**2)
        stx /= stn
        sty /= stn


        etx = e_tangent[0]
        ety = e_tangent[1]
        etn = math.sqrt(etx**2+ety**2)
        etx /= -etn
        ety /= -etn

        px = (((ey-sy)*stx+sty*sx)*etx-ety*ex*stx)/(etx*sty-ety*stx)
        py = (((-ex+sx)*sty-stx*sy)*ety+etx*ey*sty)/(etx*sty-ety*stx)

        self.control = px[0] + 1j*py[0]

    def __repr__(self):
        params = (self.start, self.control, self.end, self.weights)
        return ("RationalBezier(start={}, control={}, end={}, w={})".format(*params))

    def __eq__(self, other):
        if not isinstance(other, RationalBezier):
            return NotImplemented
        return self.start == other.start and self.end == other.end \
            and self.control == other.control \
            and self.weights == other.weights

    def __ne__(self, other):
        if not isinstance(other, RationalBezier):
            return NotImplemented
        return not self == other

    def point(self, t):
        b0=(1-t)**2
        b1=2*(1-t)*t
        b2=t**2

        c0x = self.start.real
        c0y = self.start.imag

        c1x = self.control.real
        c1y = self.control.imag

        c2x = self.end.real
        c2y = self.end.imag

        denom = b0*self.weights[0]+b1*self.weights[1]+b2*self.weights[2]

        vx = (b0*c0x*self.weights[0]+b1*c1x*self.weights[1]+b2*c2x*self.weights[2])/denom
        vy = (b0*c0y*self.weights[0]+b1*c1y*self.weights[1]+b2*c2y*self.weights[2])/denom

        return vx + 1j*vy


SAMPLE_MAX_DEPTH = 5
SAMPLE_ERROR = 1e-3


def complex_to_point(p):
    return [p.real, p.imag]


def complex_to_vect(p):
    return [p.real, p.imag]


def compute_samples(curve, start=0, end=1, error=SAMPLE_ERROR, max_depth=SAMPLE_MAX_DEPTH, depth=0):
    return compute_samples_aux(curve, start, end, curve.point(start), curve.point(end), error, max_depth, depth)


def compute_samples_aux(curve, start, end, start_point, end_point, error, max_depth, depth):
    """Recursively approximates the length by straight lines"""
    mid = (start + end)/2
    mid_point = curve.point(mid)
    length = abs(end_point - start_point)
    first_half = abs(mid_point - start_point)
    second_half = abs(end_point - mid_point)


    length2 = first_half + second_half

    res = [start, mid, end]

    if abs(length) < 1e-10:
        return []

    if depth < 2 or ((abs(length2 - length) > error) and (depth <= max_depth)):
        depth += 1
        res1 = compute_samples_aux(curve, start, mid, start_point, mid_point, error, max_depth, depth)
        res2 = compute_samples_aux(curve, mid, end, mid_point, end_point, error, max_depth, depth)
        return sorted(set(res1 + res2))
    else:
        # This is accurate enough.
        return res


def convert_svg(input_svg):
    doc = svg.Document(input_svg)
    paths = doc.flatten_all_paths()

    json_data = []
    vertices = []
    lines = []

    v_index = 1
    c_index = 0
    for ii in range(len(paths)):
        tmp = paths[ii]
        path = tmp.path

        print(str(ii+1) + "/" + str(len(paths)) + " n sub paths: " + str(len(path)))

        first = True
        for pieces in path:
            if isinstance(pieces, svg.path.Arc):
                if first:
                    if abs(pieces.start - path[-1].end) < 1e-10:
                        pieces.start = path[-1].end
                if abs(pieces.delta) > 120:
                    n_pices = math.ceil(abs(pieces.delta) / 120)
                    tmp = []
                    for p in range(n_pices):
                        t0 = float(p)/n_pices
                        t1 = float(p+1)/n_pices
                        tmp.append(RationalBezier(pieces, tstart=t0, tend=t1))
                    pieces = tmp
                else:
                    pieces = [RationalBezier(pieces)]
            else:
                pieces = [pieces]

            first = False

            for piece in pieces:
                ts = compute_samples(piece)

                if len(ts) <= 0:
                    continue


                param_ts = []

                iprev = None

                for param_t in ts:
                    param_ts.append(param_t)
                    p = piece.point(param_t)
                    xy = complex_to_point(p)

                    if not iprev:
                        istart = v_index

                    vertices += [[xy[0], xy[1]]]

                    if iprev:
                        lines += [[iprev-1, v_index-1]]

                    iprev = v_index

                    v_index += 1

                json_obj = {}
                json_obj["v_ids"] = list(range(istart-1, v_index-1))
                istart = None
                json_obj["paras"] = param_ts
                json_obj["curve_id"] = c_index
                c_index += 1

                is_set = False

                if isinstance(piece, svg.path.Line):
                    json_obj["type"] = "Line"
                    json_obj["start"] = complex_to_point(piece.start)
                    json_obj["end"] = complex_to_point(piece.end)

                    is_set = True
                elif isinstance(piece, svg.path.QuadraticBezier):
                    json_obj["type"] = "BezierCurve"
                    json_obj["degree"] = 2

                    json_obj["poles"] = [complex_to_point(piece.start), complex_to_point(piece.control), complex_to_point(piece.end)]

                    is_set = True
                elif isinstance(piece, svg.path.CubicBezier):
                    json_obj["type"] = "BezierCurve"
                    json_obj["degree"] = 3
                    json_obj["poles"] = [complex_to_point(piece.start), complex_to_point(piece.control1), complex_to_point(piece.control2), complex_to_point(piece.end)]

                    is_set = True
                elif isinstance(piece, RationalBezier):
                    json_obj["type"] = "RationalBezier"
                    json_obj["poles"] = [complex_to_point(piece.start), complex_to_point(piece.control), complex_to_point(piece.end)]
                    json_obj["weigths"] = [piece.weights[0], piece.weights[1], piece.weights[2]]

                    is_set = True

                if is_set:
                    json_data.append(json_obj)
                else:
                    print(type(piece))
                    print(piece)
                    assert(False)


    return np.array(vertices), np.array(lines), json.dumps(json_data, indent=1)



