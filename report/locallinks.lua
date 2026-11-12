-- In the PDF build, render links to repository-local files as monospace path text
-- (they are clickable on GitHub in the markdown, but meaningless inside a PDF).
function Link(el)
  local t = el.target or ""
  if t:match("^https?://") or t:match("^mailto:") or t:match("^#") then return el end
  local label = pandoc.utils.stringify(el.content)
  local path = t:gsub("^%./", ""):gsub("^%.%./", "")
  if label == path or label == "" then return pandoc.Code(path) end
  return { pandoc.Str(label), pandoc.Space(), pandoc.Str("("), pandoc.Code(path), pandoc.Str(")") }
end
