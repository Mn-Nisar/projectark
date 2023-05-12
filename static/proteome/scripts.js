ibaqSelect = document.getElementById("ibaq-select");
document.querySelector('#lable-free-bottun').addEventListener('click', function() {

$('.iba ,.name , #footer').not('#ibaq-select').hide();
ibaqSelect.style.display = 'block';

// document.querySelectorAll('body>*:not(#ibaq-select)').forEach(e => e.style.opacity = "0.3");


});

closeIbqBtn = document.getElementById("close-ibq-btn");
closeIbqBtn.addEventListener('click', function(){
  $('.iba ,.name, #footer').not('#ibaq-select').show();
// document.querySelectorAll('body>*:not(#ibaq-select').forEach(e => e.style.opacity = "1");
ibaqSelect.style.display = 'none';

});
