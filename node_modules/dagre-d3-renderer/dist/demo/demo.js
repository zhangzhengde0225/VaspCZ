const d3 = window.d3
const bodyElem = d3.select('body')
const jsElem = d3.select('#js')
const jsPanel = bodyElem.append('div').attr('id', 'jsPanel')
const cssElem = d3.select('#css')
const cssPanel = bodyElem.append('div').attr('id', 'cssPanel')

function setupPanel (panel, elem, title) {
  panel.append('h2').text(title)
  return panel.append('pre').append('code').text(elem.html().trim())
}

const jsCode = setupPanel(jsPanel, jsElem, 'JavaScript')
const cssCode = setupPanel(cssPanel, cssElem, 'CSS')

const hljsRoot = 'https://highlightjs.org/static'

bodyElem.append('link')
  .attr('rel', 'stylesheet')
  .attr('href', hljsRoot + '/demo/styles/xcode.css')
bodyElem.append('script')
  .attr('src', hljsRoot + '/highlight.site.pack.js')
  .on('load', function () {
    window.hljs.highlightBlock(jsCode.node())
    window.hljs.highlightBlock(cssCode.node())
  })
