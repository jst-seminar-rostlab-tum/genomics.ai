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
 *
 * @param categories An object containing available categories and their values as key-value pairs
 * @param setColorMode A function accepting a category to be colored
 * @param hide A function that will hide a category value in the UMAP
 * @param show A function that will set a category value in the UMAP visible if it was hidden
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
          title !== 'predictions'
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
