/* FlexAID∆S Design System — minimal interactions */
(function () {
  'use strict';

  // Install tabs
  const tabs = document.querySelectorAll('.install-tabs .tab');
  tabs.forEach(tab => {
    tab.addEventListener('click', () => {
      tabs.forEach(t => t.classList.remove('active'));
      document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));

      tab.classList.add('active');
      const panel = document.getElementById('tab-' + tab.dataset.tab);
      if (panel) panel.classList.add('active');
    });
  });

  // Simple accent color cycling on hero equation (fun detail from design system)
  const eq = document.querySelector('.equation');
  if (eq) {
    let i = 0;
    const colors = ['#25cae9', '#9571f8', '#f9b31d'];
    setInterval(() => {
      i = (i + 1) % colors.length;
      eq.style.transition = 'color 400ms ease';
      eq.style.color = colors[i];
      setTimeout(() => { eq.style.color = ''; }, 1200);
    }, 6500);
  }
})();