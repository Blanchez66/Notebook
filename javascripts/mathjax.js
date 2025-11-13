window.MathJax = {
  tex: {
     inlineMath: [ ["\$","\$"] ],
     inlineMath:[["\-\(","\)"]],
     inlineMath: [ ["\\(","\\)"] ],
     displayMath: [ ["\$\$","\$\$"] ],
     displayMath:[[ "\\[", "\\]" ]],
     processEscapes: true,
     processEnvironments: true
  },
  options: {
     ignoreHtmlClass: ".*|",
     processHtmlClass: "arithmatex"
  }
  };
document$.subscribe(() => { 
MathJax.typesetPromise()
})