Clear[f, d, b, n, N1, order, a, sol];

order = 10
lhs = (1 - d)/b - Sum[(-1)^n f[n, N1]/b^(n + 1), {n, 0, order+2}];
lhsd = lhs /. b -> Sum[a[i] d^i, {i, -1, order}];
lhsds = Series[lhsd, {d, 0, order+3}];

For [i = -1, i <= order, i++,
  c[i] = Coefficient[lhsds, d^(3 + i)] /. {f[0, N1] -> 1};
  ];

For [i = -1, i <= order, i++,
  sol[i] = Simplify[Solve[c[i] == 0, a[i]] [[1]]];
  For [j = i + 1, j <= order, j++,
    c[j] = Simplify[c[j] /. sol[i]];
    ]
  ];

Save["Documents/b_coeffs_f.txt", {order, sol}];

f[n_, N1_] := Sum[i^n, {i, 0, N1 - 1}]/N1^(n+1);
For[i = -1, i <= order, i++,
  sol[i] = Simplify[sol[i]];
  ];

Save["Documents/b_coeffs_N1.txt", {order, sol}];



