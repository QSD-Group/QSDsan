(function () {
  // Point this at the deployed Render endpoint. Overridable via a global set in conf.
  const ENDPOINT = window.QSDSAN_CHATBOT_ENDPOINT || "https://qsdsan-chatbot.onrender.com/ask";

  function el(tag, cls, html) {
    const node = document.createElement(tag);
    if (cls) node.className = cls;
    if (html !== undefined) node.innerHTML = html;
    return node;
  }

  // Minimal, safe-ish markdown: code fences, inline code, links, paragraphs.
  function renderMarkdown(md) {
    const esc = (s) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
    let html = esc(md);
    html = html.replace(/```(\w*)\n([\s\S]*?)```/g, (_m, _lang, code) => `<pre><code>${code}</code></pre>`);
    html = html.replace(/`([^`]+)`/g, "<code>$1</code>");
    html = html.replace(/\[([^\]]+)\]\((https?:[^)]+)\)/g, '<a href="$2" target="_blank" rel="noopener">$1</a>');
    html = html.replace(/\*(.+?)\*/g, "<em>$1</em>");
    html = html.split(/\n{2,}/).map((p) => (p.startsWith("<pre>") ? p : `<p>${p.replace(/\n/g, "<br>")}</p>`)).join("");
    return html;
  }

  function escAttr(s) {
    return String(s)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function citationsHtml(citations) {
    if (!citations || !citations.length) return "";
    const items = citations
      .map(
        (c) =>
          `<li>[${escAttr(c.n)}] <a href="${escAttr(c.url)}" target="_blank" rel="noopener">${escAttr(c.title)}</a> <span class="qsd-src">(${escAttr(c.source)})</span></li>`
      )
      .join("");
    return `<ul class="qsd-citations">${items}</ul>`;
  }

  function init() {
    const button = el("button", "qsd-chat-button", "Ask the docs");
    button.setAttribute("aria-label", "Open the QSDsan docs assistant");
    const panel = el("div", "qsd-chat-panel qsd-hidden");
    panel.innerHTML = `
      <div class="qsd-chat-header">QSDsan / EXPOsan docs assistant
        <button class="qsd-close" aria-label="Close">&times;</button></div>
      <div class="qsd-chat-log"></div>
      <form class="qsd-chat-form">
        <input class="qsd-chat-input" type="text" placeholder="Ask about QSDsan or EXPOsan..." autocomplete="off"/>
        <button type="submit">Send</button>
      </form>`;
    document.body.appendChild(button);
    document.body.appendChild(panel);

    const log = panel.querySelector(".qsd-chat-log");
    const form = panel.querySelector(".qsd-chat-form");
    const input = panel.querySelector(".qsd-chat-input");

    button.addEventListener("click", () => panel.classList.toggle("qsd-hidden"));
    panel.querySelector(".qsd-close").addEventListener("click", () => panel.classList.add("qsd-hidden"));

    function addMessage(cls, html) {
      const m = el("div", "qsd-msg " + cls, html);
      log.appendChild(m);
      log.scrollTop = log.scrollHeight;
      return m;
    }

    form.addEventListener("submit", async (e) => {
      e.preventDefault();
      const question = input.value.trim();
      if (!question) return;
      addMessage("qsd-user", renderMarkdown(question));
      input.value = "";
      const pending = addMessage("qsd-bot", "<em>Searching the docs...</em>");
      try {
        const resp = await fetch(ENDPOINT, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ question }),
        });
        if (!resp.ok) throw new Error("HTTP " + resp.status);
        const data = await resp.json();
        pending.innerHTML = renderMarkdown(data.answer || "") + citationsHtml(data.citations);
      } catch (err) {
        pending.innerHTML = "<em>Sorry, the docs assistant is unavailable right now.</em>";
      }
      log.scrollTop = log.scrollHeight;
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
