import { defineConfig } from '@playwright/test';
export default defineConfig({ testDir: './tests/e2e', timeout: 30000, use: { baseURL: 'http://127.0.0.1:4173', headless: true, channel: process.env.CI ? undefined : 'msedge' }, webServer: { command: 'npm run preview -- --port 4173', url: 'http://127.0.0.1:4173', reuseExistingServer: !process.env.CI } });
