module PolynomialMatrices

using DataStructures
using Polynomials
using Compat

# Import functions for overloading
import Polynomials: coeffs, degree
import Base: start, next, done
import Base: promote_rule, convert, size, length
import Base: +,-,*
import Base: getindex
import Base: copy
import Base: transpose, ctranspose
import Base: summary
import Base: insert!

# Export
export PolyMatrix
export col_degree, row_degree
export high_col_deg_matrix, high_row_deg_matrix
export is_col_proper, is_row_proper
export colred, rowred
export ltriang, rtriang

# Include files
include("polymatrix.jl")
include("arithmetic.jl")
include("reductions.jl")

end # module
