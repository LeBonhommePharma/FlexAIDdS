// main.tsx — BonhommeViewer entry point
//
// Mounts <App /> into the #root element of index.html, which is the Vite
// build entry for this app.
//
// Copyright 2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0

import { createRoot } from 'react-dom/client';
import { App } from './App';

const container = document.getElementById('root');
if (!container) {
  throw new Error('BonhommeViewer: #root element missing from index.html');
}

createRoot(container).render(<App />);
