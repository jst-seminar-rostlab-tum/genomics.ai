import { Box, ButtonBase, Stack, Typography } from "@mui/material"
import { Search as SearchIcon, FilterAlt } from "@mui/icons-material"
import { styled } from "@mui/system"
import { colors } from "shared/theme/colors"
import { useEffect, useRef, useState } from "react"

const SearchInput = styled("input")(({ theme }) => ({
  border: "none",
  outline: "none",
  fontSize: 20,
  width: "100%",
  color: theme.palette.primary.main
}))

const FilterButton = ({ onClick }) => {
  return (
    <ButtonBase
      onClick={onClick}
      disableRipple
      sx={{
        p: "12px",
        ":hover": {
          color: "primary.light"
        }
      }}
    >
      <FilterAlt />
      <Typography fontSize={20}>Filter</Typography>
    </ButtonBase>
  )
}

const Search = ({ filterComponent }) => {

  const [active, setActive] = useState(false)
  const [filterEnabled, setFilterEnabled] = useState(false)
  const filterBox = useRef()

  useEffect(() => {
    const handleFilterClose = (e) => {
      if (filterBox.current && !filterBox.current.contains(e.target)) {
        setFilterEnabled(false)
      }
    }
    window.addEventListener("click", handleFilterClose, true)
    return () => {
      window.removeEventListener("click", handleFilterClose, true);
    }
  }, [])

  return (
    <Stack direction="row" justifyContent="space-between" sx={{ position: "relative", p: "10px", border: `2px solid ${active ? colors.primary[400] : colors.primary[700]}`, borderRadius: "40px" }}>
      {/* Left part */}
      <Stack direction="row" alignItems="center" gap="5px" sx={{ marginLeft: "20px", width: "100%" }}>
        <SearchIcon sx={{ color: active ? "primary.light" : "primary.main" }} />
        <SearchInput placeholder="Search" onMouseEnter={() => setActive(true)} onMouseLeave={() => setActive(false)} />
      </Stack>
      {/* Right part */}
      <Stack direction="row" sx={{ marginRight: "20px" }}>
        <FilterButton onClick={() => setFilterEnabled(true)} />
      </Stack>
      {
        !filterEnabled ? null :
          <Box
            ref={filterBox}
            sx={{
              position: "absolute",
              zIndex: "20",
              top: "55px",
              right: "40px",
            }}
          >
            {filterComponent}
          </Box>
      }
    </Stack>
  )
}

export default Search