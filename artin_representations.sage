def Load(*s):
    for name in s:
        try:
            load(name+'.sage')
        except IOError:
            try:
                load(name+'/'+name+'.sage')
            except IOError:
                load('../'+name+'/'+name+'.sage')

Load('class_functions','artin_periods','recurrence')

def compositum(K,L):
    """
    If K and L are number fields, compositum(K,L) is a triple
    (J,f,g), where J is a number field and f:K-->J, g:L-->J
    """
    return K.composite_fields(L,both_maps=True)[0][:-1]

def restrict_element(h,f):
    # f: K --> L, h in Gal(L/Q), restrict(h,f) in Gal(K/Q)
    # is the restriction of h along f.
    K,L = f.domain(),f.codomain()
    if not K.is_galois() or not L.is_galois():
        raise ValueError("The domain and codomain of f must be Galois over Q")
    x = K.gens()[0]
    G,H = K.galois_group(),L.galois_group()
    assert h in H
    for g in G:
        if L(h(f(x))) == L(f(g(x))):
            return g
    # This should not happen
    print parent(h(f(x)))
    print parent(f(g(x)))
    print x
    print f
    print f.domain()
    print f.codomain()
    print h
    raise Exception("Could not restrict")

class ArtinReps():
    """
    Category of representations of Artin type
    For now we'll insist base is QQ
    """
    def __init__(self,coef=QQ):
        base = QQ
        if base is QQ:
            self._base = PolynomialRing(QQ,'x').one().splitting_field('x')
        else:
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
        if isinstance(s,Integer):
            self.dimension = s
            self.field = parent.base()
            self.group = self.field.galois_group()
            self.rep = {g:identity_matrix(parent.coef(),s) for g in self.G()}
        elif isinstance(s,ArtinRep):
            self.field,self.group,self.dimension = s.field,s.group,s.dimension
            self.rep = dict(s.rep)
        elif s.is_field():
            G = s.galois_group()
            self.field = s
            self.group = G
            n = G.cardinality()
            self.dimension = n
            self.rep = {}
            for g in G:
                M = [[1 if g*list(G)[j]==list(G)[i] else 0 for j in range(n)] for i in range(n)]
                self.rep[g] = Matrix(M)
            assert self.check()
            
    def G(self):
        return self.group
    
    def fill(self):
        """
        Takes the values of representation on generators and extends to the whole group
        """
        G,rep = self.group,self.rep
        while len(rep) < G.cardinality():
            b = True
            for g in G:
                for h in H:
                    a = g*h
                    if a not in rep:
                        rep[a] = rep[g] * rep[h]
                        b = False
            if b:
                raise ValueError("Set does not generate G")
        if not self.check(): # Make sure data is consistent
            raise ValueError("Values do not extend to a representation of G")
        return None
    
    def check(self):
        """
        Checks if representation is actually a group homomorphism
        """
        G = self.group
        D = self.rep
        for g in G:
            for h in G:
                if D[g*h] != D[g]*D[h]:
                    return False
        return True
    
    def character(self):
        G = self.group
        C = ClassFunctions(self.parent.coef())
        D = {}
        for g in self.rep:
            D[G.conjugacy_class(g)] = self.rep[g].trace()
        return C(D)
            
    def L_functional_equation(self):
        pass
    
    def __repr__(self):
        return str(self.rep)
        return str(self.group.as_finitely_presented_group())

    def clean(self):
        pass

    def extend_field(self,f):
        T = self.parent()
        L = f.codomain()
        T.field = L
        T.dimension = self.dimension
        T.group = L.galois_group()
        T.rep = {h:self.rep[restrict_element(h,f)] for h in T.group}
        #assert T.check()
        return T
        

    def __add__(self, s):
        C = self.parent
        T = C()
        if isinstance(s,ArtinRep):
            K,f,g = compositum(self.field,s.field)
            T.field = K
            T.group = K.galois_group()
            T.dimension = self.dimension + s.dimension
            selfE,sE = self.extend_field(f),s.extend_field(g)
            assert selfE.group is sE.group
            T.rep ={g:selfE.rep[g].block_sum(sE.rep[g]) for g in T.group}
            return T
        elif isinstance(s,Integer):
            return self + self.parent(s)
        else:
            raise ValueError("Could not add")
    
    def __radd__(self,s):
        return self + s

    def __eq__(self,other):
        pass

    def __mul__(self, s):
        C = self.parent
        T = C()
        if isinstance(s,ArtinRep):
            K,f,g = compositum(self.field,s.field)
            T.field = K
            T.group = K.galois_group()
            T.dimension = self.dimension * s.dimension
            selfE,sE = self.extend_field(f),s.extend_field(g)
            assert selfE.group is sE.group
            T.rep ={g:selfE.rep[g].tensor_product(sE.rep[g]) for g in T.group}
            return T
        elif isinstance(s,Integer):
            return self * self.parent(s)
        else:
            raise ValueError("Could not tensor")
    
    def __rmul__(self,s):
        return self * s
    
    def _pow_(self,n):
        pass

    def disp(self,string=False):
        pass

    def _latex_(self):
        pass