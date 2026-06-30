.. _tutorial_ai_assisted_development:

AI-Assisted Development
=======================

This guide covers using AI coding assistants
(`Claude Code <https://claude.com/claude-code>`_ and
`Codex <https://developers.openai.com/codex>`_) to navigate and contribute to
the QSDsan/EXPOsan codebase.

.. note::
   Sections that differ by tool use a **Claude Code / Codex** tab. The tabs are
   synced, so switching one switches them all. Everything else applies to both.

   This guide was authored and verified on Windows. The macOS/Linux equivalents
   shown (keyboard shortcuts, activation commands) follow standard conventions;
   if a step differs on your system, please
   `open an issue <https://github.com/QSD-Group/QSDsan/issues>`_.


Why AI Tools for Contributors?
------------------------------

Codebases like QSDsan combine domain-specific logic, mathematical modeling, and
software architecture across hundreds of files. AI tools help you navigate this
quickly (tracing sanunit inheritance, understanding BioSTEAM integration,
reading unfamiliar code) so you can focus your time on actual contributions
rather than orientation.

.. warning::
   **AI makes mistakes. Your review is essential.** AI tools can produce
   confident, fluent output that is wrong: a misremembered API, an outdated
   pattern, a subtle logic error. Treat everything the AI generates as a draft:
   read it, test it, and verify it against the actual code before you commit or
   open a pull request. You remain responsible for every line you contribute.


Setup
-----

.. tab-set::

   .. tab-item:: Claude Code
      :sync: cc

      .. admonition:: Before you begin

         Using Claude Code is **not free**: it requires a paid Claude plan
         (Pro or Max) or Anthropic API billing. A free account on its own will
         not work. Set this up at `claude.ai <https://claude.ai>`_ before your
         first session; account creation also requires email verification.

      #. Install VS Code from `code.visualstudio.com <https://code.visualstudio.com>`_.
         Accept all defaults; on Windows, leave all checkboxes checked.
      #. Open Extensions (``Ctrl+Shift+X`` on Windows/Linux, ``Cmd+Shift+X`` on
         Mac), search **Claude Code**, click **Install**. The button changes to
         **Disable** when done.
      #. Press ``Ctrl+L`` (Windows/Linux) or ``Cmd+L`` (Mac) to open the Claude
         Code panel and sign in.
      #. Open your workspace folder: **File → Open Folder**. See *Development
         environment* below for the recommended layout; open the folder that
         holds both repositories.

      .. note::
         This guide covers only the VS Code extension. For the terminal CLI
         installation, see the official
         `Claude Code documentation <https://code.claude.com/docs>`_ from Anthropic.

   .. tab-item:: Codex
      :sync: codex

      .. admonition:: Before you begin

         Codex access is included with eligible ChatGPT plans (such as Plus,
         Pro, Business, Edu, and Enterprise) or can be used with an OpenAI API
         key. Free or trial access may change over time, so confirm your access
         before your first session.

      #. Install VS Code from `code.visualstudio.com <https://code.visualstudio.com>`_.
         Accept all defaults; on Windows, leave all checkboxes checked.
      #. Open Extensions (``Ctrl+Shift+X`` on Windows/Linux, ``Cmd+Shift+X`` on
         Mac), search **Codex**, and install the official OpenAI Codex extension.
      #. Open your workspace folder: **File → Open Folder**. See *Development
         environment* below for the recommended layout; open the folder that
         holds both repositories.
      #. Open Codex from the VS Code sidebar. The first time, follow the prompt
         to sign in with your ChatGPT account or an API key.

      .. note::
         This guide covers only the VS Code extension. For the terminal CLI
         installation, see the official
         `Codex documentation <https://developers.openai.com/codex>`_ from OpenAI.


Development environment
-----------------------

The AI tools above help you read and change code, but you still need a working
QSDsan environment to run and test it. Follow the
:doc:`QSDsan contributing guide </CONTRIBUTING>` to fork, clone, and install
QSDsan in editable mode with its development dependencies. This is a one-time
setup.

If you also plan to contribute to ``EXPOsan``, follow the
:ref:`side-by-side workspace layout <also-contributing-to-exposan>` in the
contributing guide instead; with both packages installed editable in one shared
environment, AI tools can navigate and refactor across the two repos in a
single session, and changes in ``QSDsan`` are immediately visible to ``EXPOsan``.


Key Concepts for Contributors
-----------------------------

The AI tools, companies, and models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A quick map of the names you'll come across: the companies, their AI
assistants, the coding tools, and the models that power them.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Term
     - What it means
   * - **Anthropic & OpenAI**
     - The two AI companies whose tools this guide covers. **Anthropic** builds
       Claude; **OpenAI** builds ChatGPT. Each also makes a coding tool (below).
   * - **Claude & ChatGPT**
     - The companies' general-purpose AI assistants, the chat products you may
       already have used in a web browser. **Claude** is Anthropic's;
       **ChatGPT** is OpenAI's.
   * - **Claude Code & Codex**
     - The AI *coding* tools this guide uses. **Claude Code** (from Anthropic,
       powered by Claude) runs inside VS Code; **Codex** (from OpenAI, powered
       by its models) is available as a CLI, IDE extension, web app, and cloud
       coding agent. This guide uses the VS Code extension path. Unlike the chat
       products, these can read and edit the files in your project directly.
   * - **Model**
     - The underlying AI that does the thinking. Models come in tiers that trade
       speed against capability: Anthropic's Claude has **Haiku** (fastest,
       lightest), **Sonnet** (balanced), and **Opus** (most capable); OpenAI
       offers a similar range. Heavier models handle hard problems better but
       cost more and respond slower; both tools let you choose which to use.
   * - **Context window**
     - How much text the AI can hold in mind at once: the files it has read,
       the conversation so far, and its instructions. It is large but limited;
       on a big codebase, point the AI at the relevant files rather than
       expecting it to absorb everything, and start a fresh conversation for
       each new task.
   * - **Tokens (input & output)**
     - The unit AI text is measured in. One token is roughly three-quarters of
       a word.
       **Input tokens** are what you send the AI (your prompt, files, and the
       conversation so far); **output tokens** are what it writes back. With API
       billing they are charged separately, output usually costing more. A long
       conversation costs more because its whole history is re-sent as input
       each turn.
   * - **Knowledge cutoff**
     - Each model is trained on data up to a fixed date and knows nothing after
       it. It may be unaware of a recent library release or API change, a
       common source of confident but outdated answers. When in doubt, check
       against current documentation or point the AI at the actual code.
   * - **Agent / agentic**
     - An **agentic** tool does more than answer questions; it takes actions on
       your behalf across several steps: reading and editing files, running
       commands and tests, searching the project. Claude Code and Codex are both
       agentic, which is why they can carry out a task rather than only describe
       it.
   * - **Skill**
     - A pre-loaded instruction set that tells the AI how to approach a specific
       type of task. For example, a ``brainstorming`` skill walks it through a
       structured design process before any code is written. Both tools support
       skills; see :ref:`Skills in practice <skills-superpowers>` below.
   * - **MCP**
     - **MCP** (Model Context Protocol) is an open standard that lets AI coding
       tools connect to outside systems such as databases, documentation,
       issue trackers, and web APIs. An MCP *server* exposes one such source,
       extending what the AI can see and do beyond your project files; see
       :ref:`MCP in practice <mcp-notebooklm>` below.
   * - **Hallucination**
     - When an AI states something false with full confidence: a function that
       does not exist, a misremembered API, an invented citation. Hallucinations
       look just like correct answers, which is why you must verify the AI's
       output rather than trust it.
   * - **Subscription vs. API access**
     - Two ways to pay for the same models. A **subscription** (Claude Pro/Max,
       ChatGPT Plus/Pro/Business/Edu/Enterprise) is a flat monthly fee with a
       capped usage allowance, predictable and simplest for individuals.
       **API access** is pay-as-you-go billing through the company's developer
       platform: charged per token, with no flat fee and no hard cap. Access
       rules and free trials change over time, so check the current pricing page.
   * - **Usage / rate limits**
     - Even paid subscriptions cap how much you can use within a rolling time
       window. If you reach the limit, you wait until it resets or upgrade to a
       higher plan. Heavier models and long conversations use up the allowance
       faster.

AGENTS.md and CLAUDE.md
~~~~~~~~~~~~~~~~~~~~~~~

Text files in the repository root that give the AI project-specific
instructions: coding conventions, architecture notes, things to avoid. The
tool reads them automatically and adjusts its behavior; you don't need to do
anything, but be aware that the AI's responses are shaped by these files.

``AGENTS.md`` is the shared, tool-neutral instruction file used by Codex and
many other coding agents.

.. tab-set::

   .. tab-item:: Claude Code
      :sync: cc

      Claude Code reads ``CLAUDE.md`` if present; it may also use ``AGENTS.md``
      depending on setup, so keep shared project rules there when possible and
      Claude-only details in ``CLAUDE.md``.

   .. tab-item:: Codex
      :sync: codex

      Codex reads ``AGENTS.md``. For personal or machine-wide Codex preferences,
      use Codex configuration files rather than adding them to the repository.

In QSDsan, ``AGENTS.md`` describes the sanunit namespace structure, import
conventions, and how BioSTEAM-style units differ from QSDsan-native units.


.. _skills-superpowers:

Skills in practice: Superpowers
-------------------------------

**Superpowers** is a good example of skills in action: a plugin that bundles a
curated set of skills: brainstorming, planning, test-driven development,
debugging, code review, and more. Installing it gives you a consistent,
structured AI workflow out of the box, and it works with both Claude Code and
Codex.

Installing Superpowers
~~~~~~~~~~~~~~~~~~~~~~~

.. tab-set::

   .. tab-item:: Claude Code
      :sync: cc

      Open a Claude Code chat and run:

      .. code-block:: text

         /plugin install superpowers@claude-plugins-official

      Alternatively, install it from the web at
      `claude.com/plugins/superpowers <https://claude.com/plugins/superpowers>`_.
      Restart Claude Code if prompted.

      **Verify:** Open a new Claude Code chat, type ``/``, and confirm that
      superpowers skills appear in the command menu (look for ``brainstorming``,
      ``writing-plans``, ``test-driven-development``).

   .. tab-item:: Codex
      :sync: codex

      The Codex VS Code extension cannot install plugins itself; plugin
      installation is done from the **Codex CLI** or the **Codex App**. Because
      the extension shares its configuration with the CLI, this is a one-time
      step, and the plugin then works in the extension too.

      Install the Codex CLI (see the
      `Codex documentation <https://developers.openai.com/codex>`_), run
      ``codex``, and open the plugin browser:

      .. code-block:: text

         /plugins

      Search for ``superpowers`` and select **Install Plugin**. The Codex App
      offers the same from its **Plugins** sidebar. Restart VS Code afterward so
      the extension picks up the new plugin.

      **Verify:** in the Codex panel in VS Code, start a new session and type
      ``$``; the superpowers skills should appear in the list (look for
      ``brainstorming``, ``writing-plans``, ``test-driven-development``).

Using a skill: invoke it the way your tool expects, then describe the task:

.. tab-set::

   .. tab-item:: Claude Code
      :sync: cc

      .. code-block:: text

         /brainstorming I want to add a new SanUnit for UV disinfection.

   .. tab-item:: Codex
      :sync: codex

      .. code-block:: text

         $brainstorming I want to add a new SanUnit for UV disinfection.


.. _mcp-notebooklm:

MCP in practice: NotebookLM
---------------------------

`Google NotebookLM <https://notebooklm.google.com>`_ is a research tool where
you upload documents, papers, and notes and ask questions about them. Connecting
it through **MCP** lets Claude Code or Codex draw on your NotebookLM notebooks
while it works (for example, checking a modeling decision against the papers
behind a QSDsan unit). It is a good first MCP server to try.

.. admonition:: What you'll need

   Python 3.7+ (you almost certainly have it already for QSDsan development), a
   Google account with `NotebookLM <https://notebooklm.google.com>`_ access, and
   a web browser. The connection is made by a community tool,
   `notebooklm-mcp-cli <https://github.com/jacob-bd/notebooklm-mcp-cli>`_;
   review what it does before granting it access to your Google account.

1. Install the NotebookLM MCP tool
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a terminal, run:

.. code-block:: text

   pip install notebooklm-mcp-cli

This gives you two commands: ``nlm`` (a command-line tool) and ``notebooklm-mcp``
(the MCP server your assistant connects to).

2. Sign in to NotebookLM
~~~~~~~~~~~~~~~~~~~~~~~~~

Run ``nlm login``. This opens your browser; sign in to your Google account, and
the tool captures the sign-in automatically. Confirm it worked with
``nlm login --check``.

3. Connect it to your assistant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. tab-set::

   .. tab-item:: Claude Code
      :sync: cc

      Run:

      .. code-block:: text

         nlm setup add claude-code

   .. tab-item:: Codex
      :sync: codex

      Run:

      .. code-block:: text

         nlm setup add codex

Then check the connection with ``nlm doctor``.

4. Restart and try it
~~~~~~~~~~~~~~~~~~~~~~

.. tab-set::

   .. tab-item:: Claude Code
      :sync: cc

      Restart Claude Code, or type ``/mcp`` in the chat to reconnect.

   .. tab-item:: Codex
      :sync: codex

      Restart Codex so it picks up the new connection.

You can now ask your assistant things like:

.. code-block:: text

   List all my NotebookLM notebooks.

.. code-block:: text

   Summarize the key findings from the sources in my "QSDsan literature" notebook.

.. note::
   Setup commands can change between versions; these follow
   `notebooklm-mcp-cli v0.6.9 <https://github.com/jacob-bd/notebooklm-mcp-cli/tree/v0.6.9>`_.
   If something doesn't work, check that project's current README.


QSDsan/EXPOsan Architecture Tour with AI
----------------------------------------

Start the assistant in the QSDsan folder. Use these prompts in sequence to build
a mental model of the codebase before touching any code.

Orientation
~~~~~~~~~~~

.. code-block:: text

   Give me a high-level overview of how QSDsan is organized. What are the main modules and what does each one do?

.. code-block:: text

   What is a SanUnit? How does it relate to BioSTEAM's Unit class?

.. code-block:: text

   Show me the simplest SanUnit implementation you can find in this codebase and walk me through it.

Data flow
~~~~~~~~~

.. code-block:: text

   How does a WasteStream flow between two SanUnits in a system? Trace the data from one unit's output to the next unit's input.

.. code-block:: text

   What's the difference between a static SanUnit and a dynamic one? Show me an example of each.

Finding your way around
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   I want to understand how sedimentation is modeled. Which files should I read?

.. code-block:: text

   Where would I add a new SanUnit for UV disinfection? Walk me through the conventions I'd need to follow.


Contributor Workflow
--------------------

1. Understand the issue
~~~~~~~~~~~~~~~~~~~~~~~~

Paste the GitHub issue text into the assistant and ask:

.. code-block:: text

   What part of the codebase does this issue relate to? Which files are most likely involved?

2. Explore relevant code
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   Show me the files I'd need to read or modify to address this issue.

.. code-block:: text

   Are there existing tests that cover this behavior? Where are they?

3. Implement
~~~~~~~~~~~~

Make your change. Then ask:

.. code-block:: text

   Review this change against the conventions in AGENTS.md. Does anything look off?

4. Self-review
~~~~~~~~~~~~~~~

.. code-block:: text

   What edge cases haven't I handled in this change?

.. code-block:: text

   Is there anything in this diff that could break existing behavior?

5. Commit and open a PR
~~~~~~~~~~~~~~~~~~~~~~~~

Follow the contribution guidelines in ``QSDsan/AGENTS.md``.


Going Deeper
------------

- **QSDsan-specific skills:** The repo includes a ``qsdsan-exposan-architecture``
  skill that teaches the AI the full package structure and import conventions.
  The tool loads it when you're working on architecture-related tasks; mention
  it explicitly if you need it (``/qsdsan-exposan-architecture`` in Claude Code,
  ``$qsdsan-exposan-architecture`` in Codex).
- **Where skills live:** QSDsan currently keeps mirrored skills in
  ``QSDsan/.claude/skills/`` and ``QSDsan/.codex/skills/``. Current Codex
  documents ``.agents/skills/`` as the standard repository-scoped skill location,
  so if you want automatic Codex discovery in a fresh setup, add or symlink the
  QSDsan skills there and keep the copies synchronized.
- **Writing your own skills:** A skill is a Markdown file with structured
  instructions; no code required. See the existing skills in the directories
  above for examples.
- **Customizing AGENTS.md:** Add notes to ``AGENTS.md`` (in a fork) to tune AI
  behavior for your workflow. Changes take effect immediately on the next session.

----

*This guide was itself created with AI-assisted coding, drafted using Claude
Code, then reviewed and edited by its authors.*
