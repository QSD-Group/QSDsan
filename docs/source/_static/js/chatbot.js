(function () {
  // Point this at the deployed Render endpoint. Overridable via a global set in conf.
  const ENDPOINT = window.QSDSAN_CHATBOT_ENDPOINT || "https://qsdsan-chatbot.onrender.com/ask";

  const BOT_SVG =
    '<svg viewBox="0 0 24 24" width="20" height="20" fill="none" stroke="currentColor" ' +
    'stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">' +
    '<rect x="4" y="8" width="16" height="11" rx="2.5"/><path d="M12 8V5"/>' +
    '<circle cx="12" cy="3.4" r="1.4"/><path d="M9 13h.01M15 13h.01"/><path d="M9.5 16h5"/></svg>';
  const SEND_SVG =
    '<svg viewBox="0 0 24 24" width="18" height="18" fill="none" stroke="currentColor" ' +
    'stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">' +
    '<path d="M22 2 11 13"/><path d="M22 2 15 22l-4-9-9-4 20-7z"/></svg>';
  const COPY_SVG =
    '<svg viewBox="0 0 24 24" width="14" height="14" fill="none" stroke="currentColor" ' +
    'stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">' +
    '<rect x="9" y="9" width="11" height="11" rx="2"/><path d="M5 15V5a2 2 0 0 1 2-2h10"/></svg>';

  function el(tag, cls, html) {
    const node = document.createElement(tag);
    if (cls) node.className = cls;
    if (html !== undefined) node.innerHTML = html;
    return node;
  }

  function escapeHtml(s) {
    return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }

  function escAttr(s) {
    return escapeHtml(s).replace(/"/g, "&quot;");
  }

  // Inline formatting on already-escaped text: code, bold, italic, links.
  function renderInline(text) {
    text = text.replace(/`([^`]+)`/g, "<code>$1</code>");
    text = text.replace(/\*\*([^*]+)\*\*/g, "<strong>$1</strong>");
    text = text.replace(/(^|[^*])\*([^*\n]+)\*(?!\*)/g, "$1<em>$2</em>");
    text = text.replace(
      /\[([^\]]+)\]\((https?:[^)\s]+)\)/g,
      '<a href="$2" target="_blank" rel="noopener">$1</a>'
    );
    return text;
  }

  // Small block-level Markdown renderer: fenced code, headings, lists, paragraphs.
  // Everything is HTML-escaped up front, so only our own tags are introduced.
  function renderMarkdown(md) {
    const lines = escapeHtml(md).split("\n");
    let html = "";
    let i = 0;
    let inCode = false;
    let codeBuf = [];
    let listType = null;
    let listBuf = [];

    function flushList() {
      if (listType) {
        html += "<" + listType + ">" + listBuf.join("") + "</" + listType + ">";
        listBuf = [];
        listType = null;
      }
    }

    while (i < lines.length) {
      const line = lines[i];
      if (/^```/.test(line)) {
        if (!inCode) {
          flushList();
          inCode = true;
          codeBuf = [];
        } else {
          html += "<pre><code>" + codeBuf.join("\n") + "</code></pre>";
          inCode = false;
        }
        i++;
        continue;
      }
      if (inCode) {
        codeBuf.push(line);
        i++;
        continue;
      }
      const h = line.match(/^(#{1,6})\s+(.*)$/);
      if (h) {
        flushList();
        const lvl = Math.min(h[1].length, 6);
        html += "<h" + lvl + ' class="qsd-h">' + renderInline(h[2]) + "</h" + lvl + ">";
        i++;
        continue;
      }
      const ul = line.match(/^\s*[-*+]\s+(.*)$/);
      if (ul) {
        if (listType && listType !== "ul") flushList();
        listType = "ul";
        listBuf.push("<li>" + renderInline(ul[1]) + "</li>");
        i++;
        continue;
      }
      const ol = line.match(/^\s*\d+\.\s+(.*)$/);
      if (ol) {
        if (listType && listType !== "ol") flushList();
        listType = "ol";
        listBuf.push("<li>" + renderInline(ol[1]) + "</li>");
        i++;
        continue;
      }
      if (!line.trim()) {
        flushList();
        i++;
        continue;
      }
      flushList();
      const para = [line];
      i++;
      while (i < lines.length) {
        const nxt = lines[i];
        if (
          !nxt.trim() ||
          /^```/.test(nxt) ||
          /^#{1,6}\s/.test(nxt) ||
          /^\s*[-*+]\s+/.test(nxt) ||
          /^\s*\d+\.\s+/.test(nxt)
        ) {
          break;
        }
        para.push(nxt);
        i++;
      }
      html += "<p>" + para.map(renderInline).join("<br>") + "</p>";
    }
    flushList();
    if (inCode) html += "<pre><code>" + codeBuf.join("\n") + "</code></pre>";
    return html;
  }

  function citationsHtml(citations) {
    if (!citations || !citations.length) return "";
    const items = citations
      .map(function (c) {
        return (
          "<li>[" +
          escAttr(c.n) +
          '] <a href="' +
          escAttr(c.url) +
          '" target="_blank" rel="noopener">' +
          escAttr(c.title) +
          '</a> <span class="qsd-src">(' +
          escAttr(c.source) +
          ")</span></li>"
        );
      })
      .join("");
    return '<div class="qsd-cite-label">Sources</div><ul class="qsd-citations">' + items + "</ul>";
  }

  function init() {
    const button = el("button", "qsd-chat-button");
    button.setAttribute("aria-label", "Open the QSDsan docs assistant");
    button.innerHTML = BOT_SVG + '<span class="qsd-btn-label">Ask the docs</span>';

    const panel = el("div", "qsd-chat-panel qsd-hidden");
    panel.innerHTML =
      '<div class="qsd-chat-header"><span class="qsd-head-title">' +
      BOT_SVG +
      " QSDsan / EXPOsan assistant</span>" +
      '<button class="qsd-close" aria-label="Minimize">&minus;</button></div>' +
      '<div class="qsd-chat-log"></div>' +
      '<form class="qsd-chat-form">' +
      '<input class="qsd-chat-input" type="text" placeholder="Ask about QSDsan or EXPOsan..." autocomplete="off"/>' +
      '<button type="submit" aria-label="Send">' +
      SEND_SVG +
      "</button></form>";

    document.body.appendChild(button);
    document.body.appendChild(panel);

    const log = panel.querySelector(".qsd-chat-log");
    const form = panel.querySelector(".qsd-chat-form");
    const input = panel.querySelector(".qsd-chat-input");

    function addMessage(cls, html) {
      const m = el("div", "qsd-msg " + cls);
      m.appendChild(el("div", "qsd-msg-body", html));
      log.appendChild(m);
      log.scrollTop = log.scrollHeight;
      return m;
    }

    function attachCopy(msgEl, text) {
      if (!text) return;
      const copy = el("button", "qsd-copy", COPY_SVG + "<span>Copy</span>");
      copy.type = "button";
      copy.setAttribute("aria-label", "Copy this answer");
      copy.addEventListener("click", function () {
        if (!navigator.clipboard) return;
        navigator.clipboard.writeText(text).then(function () {
          const lab = copy.querySelector("span");
          lab.textContent = "Copied";
          setTimeout(function () {
            lab.textContent = "Copy";
          }, 1500);
        });
      });
      msgEl.appendChild(copy);
    }

    let greeted = false;
    function greet() {
      if (greeted) return;
      greeted = true;
      const intro =
        "Hi! I answer questions about **QSDsan** (API + tutorials) and **EXPOsan** " +
        "(example systems like BSM1, ADM, and HTL) from the documentation, with sources " +
        "for each answer.\n\n" +
        "Heads up: the first reply after a quiet spell can take 30 to 60 seconds while the " +
        "assistant wakes up. After that, answers come back quickly.";
      addMessage("qsd-bot qsd-intro", renderMarkdown(intro));
    }

    button.addEventListener("click", function () {
      panel.classList.toggle("qsd-hidden");
      if (!panel.classList.contains("qsd-hidden")) {
        greet();
        input.focus();
      }
    });
    panel.querySelector(".qsd-close").addEventListener("click", function () {
      panel.classList.add("qsd-hidden");
    });

    form.addEventListener("submit", async function (e) {
      e.preventDefault();
      const question = input.value.trim();
      if (!question) return;
      addMessage("qsd-user", escapeHtml(question));
      input.value = "";
      const pending = addMessage("qsd-bot", '<span class="qsd-typing">Searching the docs…</span>');
      try {
        const resp = await fetch(ENDPOINT, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ question: question }),
        });
        if (!resp.ok) throw new Error("HTTP " + resp.status);
        const data = await resp.json();
        const answer = data.answer || "";
        pending.querySelector(".qsd-msg-body").innerHTML =
          renderMarkdown(answer) + citationsHtml(data.citations);
        attachCopy(pending, answer);
      } catch (err) {
        pending.querySelector(".qsd-msg-body").innerHTML =
          '<span class="qsd-error">Sorry, the docs assistant is unavailable right now. ' +
          "If it is just waking up, wait a moment and try again.</span>";
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
