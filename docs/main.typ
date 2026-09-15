#import "@local/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "", header: "")

#set par(justify: true, leading: 0.65em)
#set heading(numbering: "1.1")
// #set math.equation(numbering: "(1)")

#align(center)[
  #text(size: 20pt, weight: "bold")[MizuRoute.jl]
  #v(0.4em)
  #text(size: 13pt)[面向 Julia 水文模型耦合的 mizuRoute 核心算法重实现]
  #v(0.6em)
  #text(size: 9pt)[技术说明与 API 文档 · v0.1.0]
]

#v(1em)

#outline(title: [目录], depth: 3)
#pagebreak()

#include "overview.typ"
#include "algorithms.typ"
#include "api.typ"
#include "validation.typ"

= 参考文献

#bibliography("refs.bib")
