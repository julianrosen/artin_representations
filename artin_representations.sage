def Load(*s):
    for name in s:
        try:
            load(name+'.sage')
        except IOError:
            try:
                load(name+'/'+name+'.sage')
            except IOError:
                load('../'+name+'/'+name+'.sage')

load('class_functions','artin_periods','recurrence')

def compositum(K,L):
    """
    If K and L are number fields, compositum(K,L) is a triple
    (J,f,g), where J is a number field and f:K-->J, g:L-->J
    """
    return K.composite_fields(L,both_maps=True)[0][:-1]

def restrict(C1,f):
    """
    Returns the restriction of the conjugacy class C1 along f
    """
    K,L = f.domain(),f.codomain()
    if not K.is_galois() or not L.is_galois():
        print K
        print L
        raise ValueError("The domain and codomain of f must be Galois over Q")
    x = K.gens()[0]
    G,H = K.galois_group(),L.galois_group()
    for C2 in G.conjugacy_classes():
        g = C2.an_element()
        for h in C1:
            if f(g(x)) == h(f(x)):
                return C2
    # This should not happen
    raise Exception("Could not restrict")

class ArtinReps():
    """
    Category of Artin type representations
    """
    def __init__(self,base=QQ,coef=QQ):
        self._base = base
        self._coef = coef
    
    def coef(self):
        return self._coef
    
    def base(self):
        return self._base
    
    def __repr__(self):
        return "The category of Artin type Galois representations of "+str(self.base())+" with coefficients in "+str(self.coef())
    
    def _latex_(self):
        pass
    
    def __call__(self,s=None):
        return ArtinRep(self,s)
    
    def __contains__(self,s):
        if isinstance(s,ArtinRep) and s.parent is self:
            return True
        else:
            return False
    
class ArtinRep():
    """
    Objects of this class are Galois representations of Artin type
    """
    def __init__(self,parent,s=None):
        self.parent = parent
        if s is None:
            s = 1
        is isinstance(s,Integer):
            self.dimension = s
            self.field = parent.base()
            self.group = self.field().galois_group()
            self.V = VectorSpace(QQ,3)
            self.rep = {g:self.V.endomorphism_ring().identity() for g in self.G()}
            
    def _repr_(self):
        return "An Artin type representation"

    def clean(self):
        pass

    def extend_field(self,f):
        pass

    def _add_(self, s):
        pass

    def __eq__(self,other):
        pass

    def _mul_(self, s):
        pass

    def _pow_(self,n):
        pass

    def disp(self,string=False):
        pass

    def _latex_(self):
        pass