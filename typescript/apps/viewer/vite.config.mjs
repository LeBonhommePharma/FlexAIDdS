// SPDX-License-Identifier: Apache-2.0

import { fileURLToPath } from 'node:url';

// Production build config for BonhommeViewer (`vite build`, entry index.html).
// The alias mirrors vitest.config.mjs and the tsconfig.json `paths` entry:
// @bonhomme/shared is bundled from source, so the build does not depend on
// packages/shared/dist having been compiled first.
export default {
  root: fileURLToPath(new URL('.', import.meta.url)),
  resolve: {
    alias: {
      '@bonhomme/shared': fileURLToPath(
        new URL('../../packages/shared/src/index.ts', import.meta.url),
      ),
    },
  },
};
