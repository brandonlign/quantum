"""Necessary suffix projection rank from a minimum stabilizer support."""
def griesmer_dimension(length,distance):
    if distance<1:raise ValueError('positive binary distance required')
    total=0;dimension=0
    while True:
        add=(distance+(1<<dimension)-1)//(1<<dimension)
        if total+add>length:return dimension
        total+=add;dimension+=1


def projection_rank_bound(m,r,t):
    if not all(isinstance(x,int) for x in [m,r,t]) or m<1 or r not in [1,2] or t<0:
        raise ValueError('require positive support length, local rank1/2, and nonnegative subset size')
    lower=0;details={}
    if t<m and m>4*t:
        # Avoid floating-point logarithms at the bound's integer boundary.
        u=0
        while (1<<(r+u+1))*(m-4*t)<=4*m-(1<<(r+2))*t:u+=1
        details['plotkin_max_kernel_dimension']=u
        lower=max(lower,2*t-u)
    if r==1 and m-t>m//2:
        delta=2*((m-2*t+1)//2)
        u=griesmer_dimension(m,delta)
        details.update(binary_distance=delta,binary_max_kernel_dimension=u)
        lower=max(lower,2*t-u)
    return dict(required_projection_rank=max(0,lower),**details)
