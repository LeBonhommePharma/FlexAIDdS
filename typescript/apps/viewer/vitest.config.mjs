// SPDX-License-Identifier: Apache-2.0

import { fileURLToPath } from 'node:url';

export default {
  root: fileURLToPath(new URL('.', import.meta.url)),
  resolve: {
    alias: {
      '@bonhomme/shared': fileURLToPath(
        new URL('../../packages/shared/src/index.ts', import.meta.url),
      ),
    },
  },
  test: {
    environment: 'node',
    include: ['src/__tests__/**/*.test.ts'],
  },
};
