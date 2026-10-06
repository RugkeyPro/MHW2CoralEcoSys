FROM node:22-alpine AS build
WORKDIR /app
COPY package.json package-lock.json ./
RUN npm ci
COPY index.html tsconfig.json vite.config.ts ./
COPY src ./src
COPY public ./public
RUN npm run build

FROM python:3.12-slim
WORKDIR /app
COPY --from=build /app/dist /app/demo
COPY serve_demo.py ./
EXPOSE 8080
CMD ["python", "serve_demo.py", "--host", "0.0.0.0", "--port", "8080"]
