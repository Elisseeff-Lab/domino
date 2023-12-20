// add to extra.js
function addLang( jQuery ) {
    $("div.sourceCode").each(function(i, v){
        var lang = $(this).children("pre").attr("class").split(' ').pop()
        var Lang = lang[0].toUpperCase() + lang.slice(1)
        $(this).before('<div class="codelabel ' + lang + '">' + Lang + ' code</div>' +
                        '<div class="codelabelspacer"></div>')
    })
}
$( document ).ready(addLang)