// R runs locally in a web worker; only public project files are mounted.
const cells = [...document.querySelectorAll('[data-r-lab]')];
if (cells.length) {
  let runtime, startup, busy = false, generation = 0;
  const mounted = new Set();
  const directories = new Set(['/home/web_user']);
  const config = fetch(new URL('assets/r-labs.json', document.baseURI)).then(r => {
    if (!r.ok) throw new Error(`Cannot load R steps (${r.status}).`);
    return r.json();
  });
  const banner = document.createElement('aside');
  banner.className = 'r-session-note';
  const message = document.createElement('span');
  message.textContent = 'Run and edit the R steps here. R starts on the first run; nothing runs automatically.';
  const restart = document.createElement('button');
  restart.type = 'button'; restart.textContent = 'Restart R';
  banner.append(message, restart); cells[0].before(banner);
  const lock = value => { busy = value; cells.forEach(c => c.querySelector('.r-run').disabled = value); };
  restart.addEventListener('click', async () => {
    generation++;
    const old = runtime; runtime = undefined; startup = undefined; mounted.clear();
    directories.clear(); directories.add('/home/web_user');
    if (old) await old.close();
    lock(false); message.textContent = 'R restarted. Run any step to start a fresh session.';
    cells.forEach(c => c.querySelector('.r-status').textContent = 'Ready to run');
  });
  async function start() {
    if (!startup) startup = (async () => {
      const { WebR } = await import('https://webr.r-wasm.org/v0.6.0/webr.mjs');
      runtime = new WebR({channelType: 3, interactive: false});
      await runtime.init();
      const version = await runtime.evalRString('R.version.string');
      message.textContent = `${version} · webR 0.6.0 · runs in your browser. Publication fitting versions are recorded in the workflow.`;
      return runtime;
    })().catch(error => { startup = undefined; throw error; });
    return startup;
  }
  async function mount(r, paths) {
    for (const path of paths) {
      if (mounted.has(path)) continue;
      if (!/^(data_derived|tables|cache)\/[a-zA-Z0-9_./-]+$/.test(path) || path.includes('..'))
        throw new Error('Unsupported input path.');
      const response = await fetch(new URL(path, document.baseURI));
      if (!response.ok) throw new Error(`Input missing: ${path} (${response.status}).`);
      const bytes = new Uint8Array(await response.arrayBuffer());
      let directory = '/home/web_user';
      for (const part of path.split('/').slice(0,-1)) {
        directory += '/' + part;
        if (!directories.has(directory)) {
          await r.FS.mkdir(directory); directories.add(directory);
        }
      }
      await r.FS.writeFile('/home/web_user/' + path, bytes);
      mounted.add(path);
    }
  }
  for (const cell of cells) {
    const code = cell.querySelector('textarea');
    const initial = code.value;
    const status = cell.querySelector('.r-status');
    const output = cell.querySelector('.r-output');
    const plots = cell.querySelector('.r-plots');
    const label = cell.querySelector('.r-output-label');
    cell.querySelector('.r-reset').addEventListener('click', () => {
      code.value = initial; output.textContent = ''; plots.replaceChildren();
      status.textContent = 'Ready to run'; label.textContent = 'Output appears here after running.';
    });
    cell.querySelector('.r-download').addEventListener('click', () => {
      const url = URL.createObjectURL(new Blob([code.value], {type:'text/plain;charset=utf-8'}));
      const a = document.createElement('a'); a.href = url; a.download = cell.dataset.rLab + '.R';
      a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
    });
    cell.querySelector('.r-run').addEventListener('click', async () => {
      if (busy) return;
      lock(true); const runGeneration = generation;
      output.textContent = ''; plots.replaceChildren();
      label.textContent = 'Output from this run';
      let shelter;
      try {
        status.textContent = 'Starting R — the first run downloads the runtime…';
        const r = await start();
        const definitions = await config;
        status.textContent = 'Loading inputs…';
        await mount(r, definitions[cell.dataset.rLab].assets);
        status.textContent = 'Running…';
        shelter = await new r.Shelter();
        const env = await shelter.evalR('new.env(parent = globalenv())');
        const started = performance.now();
        const result = await shelter.captureR(code.value, {env, withAutoprint: true,
          throwJsException: false, captureGraphics: {width: 1000, height: 760, pointsize: 12}});
        if (runGeneration !== generation) return;
        output.textContent = result.output.map(x => typeof x.data === 'string' ? x.data :
          (x.data?.message ?? JSON.stringify(x.data))).join('\n') || 'Completed without text output.';
        for (const bitmap of result.images) {
          const canvas = document.createElement('canvas');
          canvas.width = bitmap.width; canvas.height = bitmap.height;
          canvas.setAttribute('role', 'img'); canvas.setAttribute('aria-label', 'Plot generated by this R run');
          canvas.getContext('2d').drawImage(bitmap, 0, 0); bitmap.close(); plots.append(canvas);
        }
        const failed = result.output.some(x => x.type === 'error');
        status.textContent = failed ? 'R error — edit the code and run again.' :
          `Completed in ${((performance.now()-started)/1000).toFixed(2)} s`;
      } catch (error) {
        if (runGeneration !== generation) return;
        status.textContent = 'Could not complete this run';
        output.textContent = `${error.message}\nYou can download the code and run it locally. Use Restart R to retry.`;
      } finally {
        if (runGeneration === generation) { if (shelter) await shelter.purge().catch(() => {}); lock(false); }
      }
    });
  }
}
