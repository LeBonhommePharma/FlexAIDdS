import { chromium, webkit } from 'playwright';

const URL = process.argv[2] || 'http://127.0.0.1:8765/';
const engines = process.argv.includes('--safari-only')
  ? [{ name: 'webkit', launcher: webkit }]
  : [
      { name: 'webkit', launcher: webkit },
      { name: 'chromium', launcher: chromium },
    ];

const timeout = 60000;
let failed = false;

for (const engine of engines) {
  const browser = await engine.launcher.launch({ headless: true });
  const context = await browser.newContext({
    viewport: engine.name === 'webkit' ? { width: 390, height: 844 } : { width: 1280, height: 800 },
    isMobile: engine.name === 'webkit',
    hasTouch: engine.name === 'webkit',
    userAgent: engine.name === 'webkit'
      ? 'Mozilla/5.0 (iPhone; CPU iPhone OS 17_0 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.0 Mobile/15E148 Safari/604.1'
      : undefined,
  });
  const page = await context.newPage();

  const logs = [];
  page.on('console', (msg) => logs.push(`[${msg.type()}] ${msg.text()}`));
  page.on('pageerror', (err) => logs.push(`[pageerror] ${err.message}`));

  try {
    await page.goto(URL, { waitUntil: 'networkidle', timeout });
    await page.waitForTimeout(10000);

    const result = await page.evaluate(() => ({
      hasMolstar: !!window.molstar,
      build: window.__FLEXAID_MOLSTAR_BUILD__ || null,
      viewerExists: !!document.getElementById('molstar-viewer'),
      unavailable: document.getElementById('molstar-viewer')?.classList.contains('molstar-unavailable'),
      ready: document.getElementById('molstar-viewer')?.classList.contains('molstar-ready'),
      canvasCount: document.getElementById('molstar-viewer')?.querySelectorAll('canvas').length || 0,
      canvasSize: (() => {
        const canvas = document.getElementById('molstar-viewer')?.querySelector('canvas');
        return canvas ? { w: canvas.width, h: canvas.height, clientW: canvas.clientWidth, clientH: canvas.clientHeight } : null;
      })(),
      drugLabel: document.getElementById('drug-of-day-label')?.textContent?.trim() || '',
    }));

    const ok = result.canvasCount > 0 && result.canvasSize?.w > 0 && result.canvasSize?.h > 0 && !result.unavailable;
    console.log(JSON.stringify({ engine: engine.name, url: URL, ok, result, logs: logs.filter((l) => l.includes('molstar')).slice(-12) }, null, 2));
    if (!ok) failed = true;
  } catch (e) {
    failed = true;
    console.error(`TEST FAIL (${engine.name}):`, e.message);
    console.log('logs:', logs.slice(-15));
  } finally {
    await browser.close();
  }
}

process.exit(failed ? 1 : 0);