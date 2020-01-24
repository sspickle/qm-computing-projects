from visual import color

def SetArrowFromCN( cn, a):
    """
    SetArrowWithCN takes a complex number  cn  and an arrow object  a .
    This version assumes 'a' is a cylinder and sets the height of
    the cylinder based on the real part, and the color/radius based
    on the imaginary part. The radius is never set to less than 5% of the 
    magnitude of the complex number.

    """
    a.axis.z = cn.real
    a.radius = max(0.05*abs(cn), abs(cn.imag)/6.0)
    if cn.imag<0:
        a.color=color.blue
    else:
        a.color=color.red
