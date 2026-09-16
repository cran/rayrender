# Rebuild from the approved Blender OBJ exports; no OBJ is needed at runtime.
# Rscript tools/build-pig-assets.R /path/to/model_pig [output-rds]
args = commandArgs(trailingOnly = TRUE)
source_dir = if (length(args)) args[1] else "../model_pig"
destination = if (length(args) > 1) args[2] else "inst/extdata/pig.rds"
pool = new.env(parent = emptyenv())
pool$v = pool$n = pool$f = pool$materials = list()
intern = function(type, value) {
  items = pool[[type]]
  hit = which(vapply(items, identical, logical(1), value))
  if (length(hit)) {
    return(hit[1])
  }
  pool[[type]][[length(items) + 1L]] = value
  length(items) + 1L
}
quantize = function(x, scale) {
  matrix(as.integer(round(x * scale)), nrow = nrow(x), ncol = ncol(x))
}
read_obj = function(filename, offset = c(0, 0, 0)) {
  lines = readLines(filename)
  v = do.call(
    rbind,
    strsplit(sub("^v +", "", lines[grepl("^v ", lines)]), " +")
  )
  n = do.call(
    rbind,
    strsplit(sub("^vn +", "", lines[grepl("^vn ", lines)]), " +")
  )
  storage.mode(v) = storage.mode(n) = "double"
  v = sweep(v, 2, offset, "+")
  mtl = list()
  key = NULL
  for (line in readLines(sub("\\.obj$", ".mtl", filename))) {
    fields = strsplit(line, " +")[[1]]
    if (fields[1] == "newmtl") {
      key = fields[2]
      mtl[[key]] = list()
    } else if (fields[1] %in% c("Kd", "Ns", "d")) {
      mtl[[key]][[fields[1]]] = as.numeric(fields[-1])
    }
  }
  material_ids = vapply(mtl, function(m) intern("materials", m), integer(1))
  objects = list()
  key = NULL
  for (line in lines) {
    if (startsWith(line, "o ")) {
      key = substring(line, 3)
      objects[[key]] = list(tokens = list(), materials = integer())
    } else if (startsWith(line, "usemtl ")) {
      material = material_ids[[substring(line, 8)]]
    } else if (startsWith(line, "f ")) {
      ob = objects[[key]]
      ob$tokens[[length(ob$tokens) + 1L]] = strsplit(
        substring(line, 3),
        " +"
      )[[1]]
      ob$materials = c(ob$materials, material)
      objects[[key]] = ob
    }
  }
  lapply(names(objects), function(name) {
    ob = objects[[name]]
    tokens = unlist(ob$tokens)
    unique_tokens = unique(tokens)
    vi = as.integer(sub("/.*", "", unique_tokens))
    ni = as.integer(sub(".*//", "", unique_tokens))
    list(
      name = name,
      v = intern("v", quantize(v[vi, , drop = FALSE], 10000)),
      n = intern("n", quantize(n[ni, , drop = FALSE], 1000)),
      f = intern(
        "f",
        matrix(match(tokens, unique_tokens), ncol = 3, byrow = TRUE)
      ),
      materials = ob$materials
    )
  })
}
optimized = file.path(source_dir, "output/optimized/obj")
pivot = c(1.43, 2.62, 0)
offset = c(0, -.6, 0)
body = read_obj(file.path(optimized, "pig_body.obj"), offset)
heads = setNames(
  lapply(c("cheerful", "skeptical", "surprised", "excited"), function(x) {
    read_obj(file.path(optimized, paste0("head_", x, ".obj")))
  }),
  c("neutral", "skeptical", "surprised", "excited")
)
ski = read_obj(file.path(optimized, "pig_ski.obj"), offset)
spider = read_obj(
  file.path(source_dir, "output/spider/pig_spider.obj"),
  offset
)
select = function(parts, pattern) {
  Filter(function(p) grepl(pattern, p$name), parts)
}
# Fit the rigid head pose once, then reuse all expression heads in the ski outfit.
a = select(heads$excited, "^Snout_continuous_rounded_muzzle")[[1]]
b = select(ski, "^Snout_continuous_rounded_muzzle")[[1]]
a = pool$v[[a$v]] / 10000
b = pool$v[[b$v]] / 10000
stopifnot(identical(dim(a), dim(b)))
ac = sweep(a, 2, colMeans(a))
bc = sweep(b, 2, colMeans(b))
sv = svd(crossprod(ac, bc))
rotation = sv$u %*% t(sv$v)
translation = colMeans(b) - drop(colMeans(a) %*% rotation)
stopifnot(max(abs(sweep(a %*% rotation, 2, translation, "+") - b)) < .0003)
ski_head_transform = diag(4)
ski_head_transform[1:3, 1:3] = t(rotation)
ski_head_transform[1:3, 4] = translation
spider_body = select(spider, "^Spider_(skin|Tail_)")
spider_head = Filter(
  function(p) !grepl("^Spider_(skin|Tail_)", p$name),
  spider
)
# Spider head coordinates become local like the other independent expressions.
for (i in seq_along(spider_head)) {
  v = pool$v[[spider_head[[i]]$v]] / 10000
  spider_head[[i]]$v = intern(
    "v",
    quantize(sweep(v, 2, pivot + offset), 10000)
  )
}
assets = list(
  version = 1L,
  pivot = pivot + offset,
  floor = -.6,
  body = body,
  heads = heads,
  ski_body = select(ski, "^(Tail|Hoof|Body)_"),
  ski_gear = select(ski, "^(Ski|Scarf|Goggles|Binding|Balaclava)_"),
  ski_head_transform = ski_head_transform,
  spider_body = spider_body,
  spider_head = spider_head
)
# Drop intermediate geometry (e.g. the baked ski head) after extracting its pose.
parts = c(
  assets$body,
  unlist(assets$heads, recursive = FALSE),
  assets$ski_body,
  assets$ski_gear,
  assets$spider_body,
  assets$spider_head
)
for (type in c("v", "n", "f", "materials")) {
  used = sort(unique(unlist(lapply(parts, `[[`, type))))
  assets[[type]] = pool[[type]][used]
  remap = function(p) {
    p[[type]] = match(p[[type]], used)
    p
  }
  for (name in c(
    "body",
    "ski_body",
    "ski_gear",
    "spider_body",
    "spider_head"
  )) {
    assets[[name]] = lapply(assets[[name]], remap)
  }
  assets$heads = lapply(assets$heads, function(h) lapply(h, remap))
  parts = c(
    assets$body,
    unlist(assets$heads, recursive = FALSE),
    assets$ski_body,
    assets$ski_gear,
    assets$spider_body,
    assets$spider_head
  )
}
dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
saveRDS(assets, destination, compress = "xz", version = 2)
cat("Pig assets:", file.info(destination)$size, "bytes\n")
