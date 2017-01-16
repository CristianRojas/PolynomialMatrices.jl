# Computes the degree of each column of a polynomial matrix
function col_degree{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  max_deg = degree(p)
  num_col = size(p,2)

  k = fill(-1,num_col)
  for i = max_deg:-1:0
    v = coeffs(p)[i]
    for j = 1:num_col
      if k[j] < 0 && !all(v[:,j] .== zero(T))
        k[j] = i
      end
    end
  end
  k
end

# Computes the degree of each row of a polynomial matrix
function row_degree{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  max_deg = degree(p)
  num_row = size(p,1)

  k = fill(-1,num_row)
  for i = max_deg:-1:0
    v = coeffs(p)[i]
    for j = 1:num_row
      if k[j] < 0 && !all(v[j,:] .== zero(T))
        k[j] = i
      end
    end
  end
  k
end

# Computes the degree of each column, and the highest-column-degree
# coefficient matrix of a polynomial matrix
function high_col_deg_matrix{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  max_deg = degree(p)
  num_col = size(p,2)

  k   = fill(-1,1,num_col)
  Phc = zeros(T,size(p))
  for i = max_deg:-1:0
    c = coeffs(p)[i]
    for j = 1:num_col
      if k[j] < 0 && !all(c[:,j] .== zero(T))
        k[j] = i
        Phc[:,j] = c[:,j]
      end
    end
  end
  return k, Phc
end

# Computes the degree of each row, and the highest-row-degree
# coefficient matrix of a polynomial matrix
function high_row_deg_matrix{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  max_deg = degree(p)
  num_row = size(p,1)

  k   = fill(-1,num_row,1)
  Phr = zeros(size(p))
  for i = max_deg:-1:0
    c = coeffs(p)[i]
    for j = 1:num_row
      if k[j] < 0 && !all(c[j,:] .== zero(T))
        k[j] = i
        Phr[j,:] = c[j,:]
      end
    end
  end
  return k, Phr
end

# Determines if a polynomial matrix is column proper (or "column reduced")
is_col_proper{T,M,O}(p::PolyMatrix{T,M,O,1}) = true
function is_col_proper{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  Phc = high_col_deg_matrix(p)[2]
  return rank(Phc) == size(p,2)
end

# Determines if a polynomial matrix is row proper (or "row reduced")
function is_row_proper(p::PolyMatrix)
  Phr = high_row_deg_matrix(p)[2]
  return rank(Phr) == size(p,1)
end

# Computes the column reduced form of a polynomial matrix
# (via Wolovich's method)
# NOTE: Wolovich's method is known to be numerically unstable (Geurts and Praagman, 1996);
# It would be preferable to implement Geurts-Praagman's method
# NOTE: Should the procedure end with an error if p is not full rank, or simply
# indicate this as an output argument (a la SLICOT)?
function colred{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  N < 2 || size(p,1) ≥ size(p,2) ||
    error("colred: Polynomial matrix is not full column rank")

  # NOTE: This type conversion should be done in compile-time
  # (another overkill use for @generated ?)
  T1      = promote_type(T,Float16)
  p_temp  = T1 == T ? copy(p) : convert(PolyMatrix{T1,AbstractArray{T1,N},O,N}, p)

  c       = p_temp.coeffs          # Dictionary of coefficient matrices of p
  num_col = N < 2 ? 1 : size(p,2)  # Number of columns of p
  U       = PolyMatrix(eye(T,num_col),p.var)

  indN    = zeros(Int,num_col)  # Collection of non-zero entries of n
  while true
    k, Phc = high_col_deg_matrix(p_temp)
    nPhc   = nullspace(Phc)
    if size(nPhc,2) == zero(T)
      return p_temp, U
    end
    n = view(nPhc,:,1)  # One vector from the nullspace of Phc

    # Find all nonzero entries of n, and among those,
    # the one associated with the column of highest degree
    num_nz = 0               # Number of elements in indN
    Nmax   = 0               # Index of the column of p with highest degree
    max_temp = -1            # Maximum degree found so far
    for i = 1:num_col
      if n[i] != 0
        num_nz += 1
        indN[num_nz] = i
        if k[i] > max_temp
          Nmax = i
          max_temp = k[i]
        end
      end
    end

    # If there are less than 2 nonzero entries in indN,
    # the polynomial matrix is not full rank (because it has a column of zeros)
    if num_nz < 2
      error("colred: Polynomial matrix is not full column rank")
    end

    # Unimodular matrix Utemp
    Utemp = SortedDict(Dict{Int,AbstractMatrix{}}())
    insert!(Utemp, 0, eye(T,num_col))
    for i = 1:max_temp-minimum(k[indN[1:num_nz]])
      insert!(Utemp, i, zeros(T,num_col,num_col))
    end

    # Perform column reduction
    for j = 1:num_nz
      if j == Nmax
        continue
      end
      col = indN[j]

      # Update coefficient matrices
      for i = 0:k[col]
        c[max_temp-k[col]+i][:,Nmax] += n[col] / n[Nmax] * c[i][:,col]
      end

      # Update Utemp
      Utemp[max_temp-k[col]][col,Nmax] = n[col] / n[Nmax]
    end

    # Update unimodular transformation matrix U
    U = U*PolyMatrix(Utemp,(num_col,num_col),p.var)

    # Reset collection indN
    fill!(indN, 0)
  end
end

# Computes the simultaneous column reduced form of two polynomial matrices `p1` and `p2`,
# according to `p1`(via Wolovich's method)
function colred{T,M1,M2,O1,O2,N1,N2}(p1::PolyMatrix{T,M1,O1,N1},
  p2::PolyMatrix{T,M2,O2,N2})
  size(p1,2) == size(p2,2) || (N1 < 2 && N2 < 2) ||
    error("colred: Both polynomial matrices should have the same number of columns")
  p1.var == p2.var ||
    error("colred: Both polynomial matrices should be in the same variable")
  N1 < 2 || size(p1,1) ≥ size(p1,2) ||
    error("colred: Polynomial matrix `p1` is not full column rank")

  T1       = promote_type(T,Float16)
  p1_temp  = T1 == T ? copy(p1) : convert(PolyMatrix{T1,AbstractArray{T1,N1},O1,N1}, p1)
  p2_temp  = T1 == T ? copy(p2) : convert(PolyMatrix{T1,AbstractArray{T1,N2},O2,N2}, p2)
  c1       = coeffs(p1_temp)          # Dictionary of coefficient matrices of p1
  c2       = coeffs(p2_temp)          # Dictionary of coefficient matrices of p2
  num_col  = N1 < 2 ? 1 : size(p1,2)  # Number of columns of p1 and p2

  indN    = zeros(Int,num_col)  # Collection of non-zero entries of n
  while true
    k1, Phc = high_col_deg_matrix(p1_temp)
    k2      = col_degree(p2_temp)
    nPhc    = nullspace(Phc)
    if size(nPhc,2) == zero(T)
      return p1_temp, p2_temp
    end

    n = view(nPhc,:,1)  # One vector from the nullspace of Phc

    # Find all nonzero entries of n, and among those,
    # the one associated with the column of highest degree
    num_nz   = 0             # Number of elements in indN
    Nmax     = 0             # Index of the column of p1 with highest degree
    max_temp = -1            # Maximum degree found so far
    for i = 1:num_col
      if n[i] != 0
        num_nz += 1
        indN[num_nz] = i
        if k1[i] > max_temp
          Nmax     = i
          max_temp = k1[i]
        end
      end
    end

    # If there are less than 2 nonzero entries in indN,
    # the polynomial matrix is not full rank (because it has a column of zeros)
    if num_nz < 2
      error("colred: Polynomial matrix is not full rank")
    end

    # Perform column reduction
    for j = 1:num_nz
      if j == Nmax
        continue
      end
      col = indN[j]

      # Update coefficient matrices of p1
      for i = 0:k1[col]
        c1[max_temp-k1[col]+i][:,Nmax] += n[col] / n[Nmax] * c1[i][:,col]
      end

      # Update coefficient matrices of p2
      for i = 0:k2[col]
        c2[max_temp-k1[col]+i][:,Nmax] += n[col] / n[Nmax] * c2[i][:,col]
      end
    end

    # Reset collection indN
    fill!(indN, 0)
  end
end

# Computes the row reduced form of a polynomial matrix
# (via Wolovich's method)
function rowred{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  (N < 2 && size(p,1) ≤ 1) || size(p,1) ≤ size(p,2) ||
    error("rowred: Polynomial matrix is not full row rank")

  T1      = promote_type(T,Float16)
  p_temp  = T1 == T ? copy(p) : convert(PolyMatrix{T1,AbstractArray{T1,N},O,N}, p)
  c       = coeffs(p_temp)  # Dictionary of coefficient matrices of p
  num_row = size(p,1)      # Number of rows of p
  U       = PolyMatrix(eye(T,num_row),p.var)

  indN    = zeros(Int,num_row)  # Collection of non-zero entries of n
  while true
    k, Phr = high_row_deg_matrix(p_temp)
    nPhr   = nullspace(Phr')
    if size(nPhr,2) == zero(T)
      return p_temp, U
    end

    n = view(nPhr,:,1)  # One vector from the nullspace of Phc

    # Find all nonzero entries of n, and among those,
    # the one associated with the row of highest degree
    num_nz = 0               # Number of elements in indN
    Nmax   = 0               # Index of the row of p with highest degree
    max_temp = -1            # Maximum degree found so far
    for i = 1:num_row
      if n[i] != 0
        num_nz += 1
        indN[num_nz] = i
        if k[i] > max_temp
          Nmax = i
          max_temp = k[i]
        end
      end
    end

    # If there are less than 2 nonzero entries in indN,
    # the polynomial matrix is not full rank (because it has a row of zeros)
    if num_nz < 2
      error("rowred: Polynomial matrix is not full rank")
    end

    # Unimodular matrix Utemp
    Utemp = SortedDict(Dict{Int,AbstractMatrix{}}())
    insert!(Utemp, 0, eye(T,num_row))
    for i = 1:max_temp-minimum(k[indN[1:num_nz]])
      insert!(Utemp, i, zeros(T,num_row,num_row))
    end

    # Perform row reduction
    for j = 1:num_nz
      if j == Nmax
        continue
      end
      row = indN[j]

      # Update coefficient matrices
      for i = 0:k[row]
        c[max_temp-k[row]+i][Nmax,:] += n[row] / n[Nmax] * c[i][row,:]
      end

      # Update Utemp
      Utemp[max_temp-k[row]][Nmax,row] = n[row] / n[Nmax]
    end

    # Update unimodular transformation matrix U
    U = PolyMatrix(Utemp,(num_row,num_row),p.var)*U

    # Reset collection indN
    fill!(indN, 0)
  end
end

# Computes the simultaneous row reduced form of two polynomial matrices `p1` and `p2`,
# according to `p1`(via Wolovich's method)
function rowred{T,M1,M2,O1,O2,N1,N2}(p1::PolyMatrix{T,M1,O1,N1},
  p2::PolyMatrix{T,M2,O2,N2})
  size(p1,1) == size(p2,1) ||
    error("rowred: Both polynomial matrices should have the same number of rows")
  p1.var == p2.var ||
    error("rowred: Both polynomial matrices should be in the same variable")
  (N1 < 2 && size(p1,1) ≤ 1) || size(p1,1) ≤ size(p1,2) ||
    error("rowred: Polynomial matrix `p1` is not full row rank")

  T1       = promote_type(T,Float16)
  p1_temp  = T1 == T ? copy(p1) : convert(PolyMatrix{T1,AbstractArray{T1,N1},O1,N1}, p1)
  p2_temp  = T1 == T ? copy(p2) : convert(PolyMatrix{T1,AbstractArray{T1,N2},O2,N2}, p2)
  c1       = coeffs(p_temp1)  # Dictionary of coefficient matrices of p1
  c2       = coeffs(p_temp2)  # Dictionary of coefficient matrices of p2
  num_row  = size(p1,1)       # Number of rows of p1 and p2

  indN    = zeros(Int,num_row)  # Collection of non-zero entries of n
  while true
    k1, Phr = high_row_deg_matrix(p1_temp)
    k2      = row_degree(p2_temp)
    nPhr    = nullspace(Phr')
    if size(nPhr,2) == zero(T)
      return p1_temp, p2_temp
    end

    n = view(nPhr,:,1)  # One vector from the nullspace of Phr

    # Find all nonzero entries of n, and among those,
    # the one associated with the row of highest degree
    num_nz   = 0             # Number of elements in indN
    Nmax     = 0             # Index of the row of p1 with highest degree
    max_temp = -1            # Maximum degree found so far
    for i = 1:num_row
      if n[i] != 0
        num_nz += 1
        indN[num_nz] = i
        if k1[i] > max_temp
          Nmax     = i
          max_temp = k1[i]
        end
      end
    end

    # If there are less than 2 nonzero entries in indN,
    # the polynomial matrix is not full rank (because it has a row of zeros)
    if num_nz < 2
      error("rowred: Polynomial matrix is not full rank")
    end

    # Perform row reduction
    for j = 1:num_nz
      if j == Nmax
        continue
      end
      row = indN[j]

      # Update coefficient matrices of p1
      for i = 0:k1[row]
        c1[max_temp-k1[row]+i][Nmax,:] += n[row] / n[Nmax] * c1[i][row,:]
      end

      # Update coefficient matrices of p2
      for i = 0:k2[row]
        c2[max_temp-k1[row]+i][Nmax,:] += n[row] / n[Nmax] * c2[i][row,:]
      end
    end

    # Reset collection indN
    fill!(indN, 0)
  end
end

# Computes a triangular form for a polynomial matrix,
# based on a left unimodular matrix, and it returns both the
# triangular matrix and the unimodular transformation used
# (based on the Henrion-Sebek algorithm)
function ltriang{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  n, m = size(p)
  m ≥ n || error("ltriang: p must have full row rank")

  dp = degree(p)

  # Compute upper bound on degree of unimodular matrix
  col   = sort(col_degree(p),rev=true)
  row   = row_degree(p)
  dUmax = min(sum(col[1:n]),sum(row))

  Q       = AbstractArray
  R       = AbstractArray
  dU      = Int
  shape   = AbstractArray
  shape_b = AbstractArray

  for dU = 0:dUmax

    # Build Sylvester matrix
    Rd = zeros(T,n*(dU+1),m*(dp+dU+1))
    for (k,c) in coeffs(p)
      for i = 0:dU
        Rd[i*n+1:(i+1)*n,dp-k+i+1:dp+dU+1:dp-k+i+1+(m-1)*(dp+dU+1)] = c
      end
    end

    # Triangularize Sylvester matrix
    Q, R = qr(Rd; thin=false)

    # Extract triangular shape
    shape_b = zeros(Int,n*(dU+1),2)
    for i = 1:n*(dU+1)
      for j = 1:m*(dp+dU+1)
        if !(R[i,j] ≈ 0)
          shape_b[i,1] = j
          shape_b[i,2] = ((j -1) ÷ (dp+dU+1)) + 1
          break
        end
      end
    end

    # Choose suitable rows for triangularized p
    shape      = zeros(Int,n)
    curr_group = shape_b[1,2]
    iter       = 1
    for j = 2:n*(dU+1)

      # Terminate if all n rows have been chosen
      if iter == n
        break
      end

      if shape_b[j,2] ≠ curr_group
        shape[iter] = j-1
        curr_group  = shape_b[j,2]
        iter       += 1
      end
    end

    if shape_b[n*(dU+1),2] ≠ shape_b[n*(dU+1)-1,2]
      shape[n]  = n*(dU+1)
      iter     += 1
    end

    if iter > n
      break
    end
  end

  # Return the chosen rows of Q and R
  c_triang = SortedDict(Dict{Int,AbstractArray{eltype(Q),N}}())
  for k = 0:dp+dU
    c_triang[k] = R[shape, dp+dU-k+1:dp+dU+1:dp+dU-k+1+(m-1)*(dp+dU+1)]
  end
  p_triang = PolyMatrix(c_triang, (n,m), p.var)

  c_triang = SortedDict(Dict{Int,AbstractArray{eltype(Q),N}}())
  for k = 0:dU
    c_triang[k] = Q[(dU-k)*n+1:(dU-k+1)*n,shape]'
  end
  Q_triang = PolyMatrix(c_triang, (n,n), p.var)

  p_triang, Q_triang
end

# Computes a triangular form for a polynomial matrix,
# based on a right unimodular matrix, and it returns both the
# triangular matrix and the unimodular transformation used
# NOTE: It would be preferable to avoid the transposes, but the QR
# factorization function only applies to left factorizations
function rtriang{T,M,O,N}(p::PolyMatrix{T,M,O,N})
  p_triang, Q_triang = ltriang(p')
  return p_triang', Q_triang'
end
