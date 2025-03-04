Clear [order, f, fb, d, B, i, b]
Clear [fbs, fbs, g, gb, gbs, gbs]
Clear [lhs, h, hb, hbs, rhs]
Clear [j, l, m, eq, eqs, c, sol]

order = 10;

(* f[x_] = Sum[PolyGamma[n[i]+ x], {i,1,k2}]; *)
fb = f[Sum[b[j] * d^(j+1), {j,0,order}]];
fbs = d*Series[fb, {d,0,order-1}];


g  = -k2*PolyGamma[1+d*B];
gb = g /. {B -> Sum[b[l] * d^l, {l,0,order}]};  
gbs = d*Series[gb, {d,0,order-1}];      

lhs = Normal[fbs] + Normal[gbs];

h = k1/B +PolyGamma[B] - PolyGamma[N+B];
hb = h /. {B -> Sum[b[m] * d^m, {m,0,order}]};
hbs = Series[hb, {d,0,order}];

rhs = Normal[hbs] /. {(k1/b[0] +PolyGamma[b[0]]-PolyGamma[b[0]+N]) ->0};

eq = Simplify[lhs+rhs];

For [i = 1, i <= order, i++,
  c[i] = Coefficient[eq, d^i];
  ];

For [i = 1, i <= order, i++,
  sol[i] = Simplify[Solve[c[i] == 0, b[i]] [[1]]];
(*  For [j = i + 1, j <= order, j++,
    c[j] = c[j] /. sol[i];
    ]*)
  ];

Save["Documents/k_coeffs_b.txt", {order, sol}];

For [i = 1, i <= order, i++,  
  sol[i] = Simplify[sol[i]];
  ];

For [i = 1, i <= order, i++,
  sol[i] = Simplify[sol[i]];  
  For [j = i + 1, j <= order, j++,
    sol[j] = sol[j] /. sol[i];
    ];
  ];

Save["Documents/k_coeffs_simpl.txt", {order, sol}];
 
