import {
  ExpandLess, ExpandMore, Visibility, VisibilityOff,
} from '@mui/icons-material';
import OpacityIcon from '@mui/icons-material/Opacity';
import {
  Box, Collapse, IconButton, ListItem, ListItemButton, Tooltip, Typography,
} from '@mui/material';
import React, { useEffect, useState } from 'react';
import { colors } from 'shared/theme/colors';

const activatedColor = colors.primary['400'];
const deactivatedColor = colors.primary['200'];

/**
 * One category entry in the categories list
 * @param title Name of the cateogiry
 * @param values An object containing the values and colors of a category as key-value pairs
 * @param colored True if the category is selected as the color mode
 * @param toggleColored Function that should be executed when the color button is clicked.
 * Expects no parameter
 * @param hide A function that will hide a category value in the UMAP.
 * Expects a value as parameter
 * @param show A function that will set a category value in the UMAP visible if it was hidden
 * Expects a value as parameter
 */
function Category({
  title, values, colored, toggleColored, hide, show,
}) {
  const [open, setOpen] = useState(false);

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column' }}>
      <ListItemButton onClick={() => setOpen(!open)} sx={{ p: 0 }}>
        <ListItem
          secondaryAction={(
            <IconButton edge="end" onClick={(e) => { toggleColored(); e.stopPropagation(); }}>
              <OpacityIcon sx={{ color: colored ? activatedColor : deactivatedColor }} />
            </IconButton>
          )}
          sx={{ pl: 1 }}
        >
          {open
            ? <ExpandLess sx={{ transform: 'rotate(180deg)' }} />
            : <ExpandMore sx={{ transform: 'rotate(-90deg)' }} />}
          <Typography sx={{ pr: 1 }} noWrap>{title}</Typography>

        </ListItem>
      </ListItemButton>
      <Collapse in={open}>
        <Box sx={{ display: 'flex', flexDirection: 'column' }}>
          {Object.entries(values).map(([value, color]) => (
            <Value
              title={value}
              color={color}
              key={value}
              hide={() => hide(value)}
              show={() => show(value)}
            />
          ))}
        </Box>
      </Collapse>
    </Box>
  );
}

/**
 * One value entry belonging to a category
 * @param title Name of the value
 * @param color Color of the value in the UMAP
 * @param hide A function that will hide the category value in the UMAP.
 * Expects no parameter
 * @param show A function that will set the category value in the UMAP visible if it was hidden
 * Expects no parameter
 */
function Value({
  title, color, hide, show,
}) {
  const [visible, setVisible] = useState(true);

  return (
    <Box sx={{
      display: 'flex', alignItems: 'center', pl: 3, pr: 2,
    }}
    >
      <IconButton onClick={() => {
        if (visible) {
          hide();
        } else {
          show();
        }
        setVisible(!visible);
      }}
      >
        { visible
          ? <Visibility sx={{ color }} />
          : <VisibilityOff sx={{ color: deactivatedColor }} />}
      </IconButton>
      <Tooltip title={title} placement="right">
        <Typography sx={{ flexGrow: 1 }} noWrap>{title}</Typography>
      </Tooltip>
    </Box>
  );
}

/**
 * Displays the catogries of a result in a vertical layout
 * @param categories An object containing available categories and their values as key-value pairs
 * @param setColorMode A function accepting a category to be colored
 * @param hide A function that will hide a category value in the UMAP.
 * Expects category and value as parameter
 * @param show A function that will set a category value in the UMAP visible if it was hidden
 * Expects category and value as parameter
 */
function GeneMapperCategories({
  categories, setColorMode, hide, show,
}) {
  const [coloredCategoryTitle, setColoredCategoryTitle] = useState(null);

  const handleSetColorMode = (colorMode) => {
    setColoredCategoryTitle(colorMode);
    setColorMode(colorMode);
  };

  useEffect(() => {
    if (categories) {
      const coloringModes = Object.keys(categories);

      if (coloringModes.includes('cell_type')) {
        handleSetColorMode('cell_type');
      } else if (coloringModes.includes('batch')) {
        handleSetColorMode('batch');
      } else if (coloringModes.length) {
        handleSetColorMode(coloringModes[0]);
      }
    }
  }, [categories]);

  return (
    <Box sx={{ maxWidth: '270px' }}>
      { categories
        ? Object.entries(categories).map(([title, values]) => (
          title !== 'predicted'
            ? (
              <Category
                key={title}
                title={title}
                values={values}
                colored={title === coloredCategoryTitle}
                toggleColored={() => handleSetColorMode(title)}
                hide={(value) => {
                  hide(title, value);
                }}
                show={(value) => {
                  show(title, value);
                }}
              />
            ) : null
        ))
        : null}
    </Box>

  );
}

export default GeneMapperCategories;
