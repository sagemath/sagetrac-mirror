{%- extends 'display_priority.tpl' -%}

{%- block header -%}
.. escape-backslashes
.. default-role:: math
{% endblock header %}


{% block in_prompt %}

{% endblock in_prompt %}

{% block output_prompt %}
{% endblock output_prompt %}

{% block input scoped%}
{%- if cell.source.strip() -%}
::

{{ cell.source | add_prompts(first='sage: ', cont='....: ') | indent}}
{%- endif -%}
{%- endblock input -%}

{% block error %}

::

{{ super() }}
{% endblock error %}

{% block traceback_line %}
{{ line | indent | strip_ansi }}
{% endblock traceback_line %}

{%- block execute_result -%}
{%- block data_priority scoped -%}
{{ super() }}
{% endblock %}
{% endblock execute_result %}

{%- block stream -%}
{{ output.text | indent }}
{% endblock stream %}

{% block data_svg %}

.. image:: {{ output.metadata.filenames['image/svg+xml'] | urlencode }}
{% endblock data_svg %}

{% block data_png %}

.. image:: {{ output.metadata.filenames['image/png'] | urlencode }}
{% endblock data_png %}

{% block data_jpg %}

.. image:: {{ output.metadata.filenames['image/jpeg'] | urlencode }}
{% endblock data_jpg %}

{% block data_markdown %}
{{ output.data['text/markdown'] | convert_pandoc("markdown", "rst") }}
{% endblock data_markdown %}

{% block data_latex %}

.. math::

{{ output.data['text/latex'] | strip_dollars | indent }}
{% endblock data_latex %}

{%- block data_text scoped -%}
{{ output.data['text/plain'] | indent }}
{% endblock data_text %}

{% block data_html scoped %}

.. raw:: html

{{ output.data['text/html'] | indent }}
{% endblock data_html %}

{% block markdowncell scoped %}

{{ cell.source | convert_pandoc("markdown", "rst") }}
{% endblock markdowncell %}

{%- block rawcell scoped -%}
{%- if cell.metadata.get('raw_mimetype', '').lower() in resources.get('raw_mimetypes', ['']) %}
{{cell.source}}
{% endif -%}
{%- endblock rawcell -%}

{% block headingcell scoped %}
{{ ("#" * cell.level + cell.source) | replace('\n', ' ') | convert_pandoc("markdown", "rst") }}
{% endblock headingcell %}

{% block unknowncell scoped %}
unknown type  {{cell.type}}
{% endblock unknowncell %}
